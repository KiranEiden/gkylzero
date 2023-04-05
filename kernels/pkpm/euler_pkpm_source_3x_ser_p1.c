#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void euler_pkpm_source_3x_ser_p1(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* out) 
{ 
  // qmem:             q/m*EM fields.
  // vlasov_pkpm_moms: [rho, p_parallel, q_parallel], Moments computed from kinetic equation in pkpm model.
  // euler_pkpm:       Input fluid variables.
  // out:              Output increment
  const double *rho = &vlasov_pkpm_moms[0]; 

  const double *Ex = &qmem[0]; 
  const double *Ey = &qmem[8]; 
  const double *Ez = &qmem[16]; 
  const double *Bx = &qmem[24]; 
  const double *By = &qmem[32]; 
  const double *Bz = &qmem[40]; 

  const double *rhoux = &euler_pkpm[0]; 
  const double *rhouy = &euler_pkpm[8]; 
  const double *rhouz = &euler_pkpm[16]; 

  double *outrhoux = &out[0]; 
  double *outrhouy = &out[8]; 
  double *outrhouz = &out[16]; 

  outrhoux[0] += (-0.3535533905932737*By[7]*rhouz[7])+0.3535533905932737*Bz[7]*rhouy[7]+0.3535533905932737*Ex[7]*rho[7]-0.3535533905932737*By[6]*rhouz[6]+0.3535533905932737*Bz[6]*rhouy[6]+0.3535533905932737*Ex[6]*rho[6]-0.3535533905932737*By[5]*rhouz[5]+0.3535533905932737*Bz[5]*rhouy[5]+0.3535533905932737*Ex[5]*rho[5]-0.3535533905932737*By[4]*rhouz[4]+0.3535533905932737*Bz[4]*rhouy[4]+0.3535533905932737*Ex[4]*rho[4]-0.3535533905932737*By[3]*rhouz[3]+0.3535533905932737*Bz[3]*rhouy[3]+0.3535533905932737*Ex[3]*rho[3]-0.3535533905932737*By[2]*rhouz[2]+0.3535533905932737*Bz[2]*rhouy[2]+0.3535533905932737*Ex[2]*rho[2]-0.3535533905932737*By[1]*rhouz[1]+0.3535533905932737*Bz[1]*rhouy[1]+0.3535533905932737*Ex[1]*rho[1]-0.3535533905932737*By[0]*rhouz[0]+0.3535533905932737*Bz[0]*rhouy[0]+0.3535533905932737*Ex[0]*rho[0]; 
  outrhoux[1] += (-0.3535533905932737*By[6]*rhouz[7])+0.3535533905932737*Bz[6]*rhouy[7]+0.3535533905932737*Ex[6]*rho[7]+0.3535533905932737*rho[6]*Ex[7]+0.3535533905932737*rhouy[6]*Bz[7]-0.3535533905932737*rhouz[6]*By[7]-0.3535533905932737*By[3]*rhouz[5]+0.3535533905932737*Bz[3]*rhouy[5]+0.3535533905932737*Ex[3]*rho[5]+0.3535533905932737*rho[3]*Ex[5]+0.3535533905932737*rhouy[3]*Bz[5]-0.3535533905932737*rhouz[3]*By[5]-0.3535533905932737*By[2]*rhouz[4]+0.3535533905932737*Bz[2]*rhouy[4]+0.3535533905932737*Ex[2]*rho[4]+0.3535533905932737*rho[2]*Ex[4]+0.3535533905932737*rhouy[2]*Bz[4]-0.3535533905932737*rhouz[2]*By[4]-0.3535533905932737*By[0]*rhouz[1]+0.3535533905932737*Bz[0]*rhouy[1]+0.3535533905932737*Ex[0]*rho[1]+0.3535533905932737*rho[0]*Ex[1]+0.3535533905932737*rhouy[0]*Bz[1]-0.3535533905932737*rhouz[0]*By[1]; 
  outrhoux[2] += (-0.3535533905932737*By[5]*rhouz[7])+0.3535533905932737*Bz[5]*rhouy[7]+0.3535533905932737*Ex[5]*rho[7]+0.3535533905932737*rho[5]*Ex[7]+0.3535533905932737*rhouy[5]*Bz[7]-0.3535533905932737*rhouz[5]*By[7]-0.3535533905932737*By[3]*rhouz[6]+0.3535533905932737*Bz[3]*rhouy[6]+0.3535533905932737*Ex[3]*rho[6]+0.3535533905932737*rho[3]*Ex[6]+0.3535533905932737*rhouy[3]*Bz[6]-0.3535533905932737*rhouz[3]*By[6]-0.3535533905932737*By[1]*rhouz[4]+0.3535533905932737*Bz[1]*rhouy[4]+0.3535533905932737*Ex[1]*rho[4]+0.3535533905932737*rho[1]*Ex[4]+0.3535533905932737*rhouy[1]*Bz[4]-0.3535533905932737*rhouz[1]*By[4]-0.3535533905932737*By[0]*rhouz[2]+0.3535533905932737*Bz[0]*rhouy[2]+0.3535533905932737*Ex[0]*rho[2]+0.3535533905932737*rho[0]*Ex[2]+0.3535533905932737*rhouy[0]*Bz[2]-0.3535533905932737*rhouz[0]*By[2]; 
  outrhoux[3] += (-0.3535533905932737*By[4]*rhouz[7])+0.3535533905932737*Bz[4]*rhouy[7]+0.3535533905932737*Ex[4]*rho[7]+0.3535533905932737*rho[4]*Ex[7]+0.3535533905932737*rhouy[4]*Bz[7]-0.3535533905932737*rhouz[4]*By[7]-0.3535533905932737*By[2]*rhouz[6]+0.3535533905932737*Bz[2]*rhouy[6]+0.3535533905932737*Ex[2]*rho[6]+0.3535533905932737*rho[2]*Ex[6]+0.3535533905932737*rhouy[2]*Bz[6]-0.3535533905932737*rhouz[2]*By[6]-0.3535533905932737*By[1]*rhouz[5]+0.3535533905932737*Bz[1]*rhouy[5]+0.3535533905932737*Ex[1]*rho[5]+0.3535533905932737*rho[1]*Ex[5]+0.3535533905932737*rhouy[1]*Bz[5]-0.3535533905932737*rhouz[1]*By[5]-0.3535533905932737*By[0]*rhouz[3]+0.3535533905932737*Bz[0]*rhouy[3]+0.3535533905932737*Ex[0]*rho[3]+0.3535533905932737*rho[0]*Ex[3]+0.3535533905932737*rhouy[0]*Bz[3]-0.3535533905932737*rhouz[0]*By[3]; 
  outrhoux[4] += (-0.3535533905932737*By[3]*rhouz[7])+0.3535533905932737*Bz[3]*rhouy[7]+0.3535533905932737*Ex[3]*rho[7]+0.3535533905932737*rho[3]*Ex[7]+0.3535533905932737*rhouy[3]*Bz[7]-0.3535533905932737*rhouz[3]*By[7]-0.3535533905932737*By[5]*rhouz[6]+0.3535533905932737*Bz[5]*rhouy[6]+0.3535533905932737*Ex[5]*rho[6]+0.3535533905932737*rho[5]*Ex[6]+0.3535533905932737*rhouy[5]*Bz[6]-0.3535533905932737*rhouz[5]*By[6]-0.3535533905932737*By[0]*rhouz[4]+0.3535533905932737*Bz[0]*rhouy[4]+0.3535533905932737*Ex[0]*rho[4]+0.3535533905932737*rho[0]*Ex[4]+0.3535533905932737*rhouy[0]*Bz[4]-0.3535533905932737*rhouz[0]*By[4]-0.3535533905932737*By[1]*rhouz[2]+0.3535533905932737*Bz[1]*rhouy[2]+0.3535533905932737*Ex[1]*rho[2]+0.3535533905932737*rho[1]*Ex[2]+0.3535533905932737*rhouy[1]*Bz[2]-0.3535533905932737*rhouz[1]*By[2]; 
  outrhoux[5] += (-0.3535533905932737*By[2]*rhouz[7])+0.3535533905932737*Bz[2]*rhouy[7]+0.3535533905932737*Ex[2]*rho[7]+0.3535533905932737*rho[2]*Ex[7]+0.3535533905932737*rhouy[2]*Bz[7]-0.3535533905932737*rhouz[2]*By[7]-0.3535533905932737*By[4]*rhouz[6]+0.3535533905932737*Bz[4]*rhouy[6]+0.3535533905932737*Ex[4]*rho[6]+0.3535533905932737*rho[4]*Ex[6]+0.3535533905932737*rhouy[4]*Bz[6]-0.3535533905932737*rhouz[4]*By[6]-0.3535533905932737*By[0]*rhouz[5]+0.3535533905932737*Bz[0]*rhouy[5]+0.3535533905932737*Ex[0]*rho[5]+0.3535533905932737*rho[0]*Ex[5]+0.3535533905932737*rhouy[0]*Bz[5]-0.3535533905932737*rhouz[0]*By[5]-0.3535533905932737*By[1]*rhouz[3]+0.3535533905932737*Bz[1]*rhouy[3]+0.3535533905932737*Ex[1]*rho[3]+0.3535533905932737*rho[1]*Ex[3]+0.3535533905932737*rhouy[1]*Bz[3]-0.3535533905932737*rhouz[1]*By[3]; 
  outrhoux[6] += (-0.3535533905932737*By[1]*rhouz[7])+0.3535533905932737*Bz[1]*rhouy[7]+0.3535533905932737*Ex[1]*rho[7]+0.3535533905932737*rho[1]*Ex[7]+0.3535533905932737*rhouy[1]*Bz[7]-0.3535533905932737*rhouz[1]*By[7]-0.3535533905932737*By[0]*rhouz[6]+0.3535533905932737*Bz[0]*rhouy[6]+0.3535533905932737*Ex[0]*rho[6]+0.3535533905932737*rho[0]*Ex[6]+0.3535533905932737*rhouy[0]*Bz[6]-0.3535533905932737*rhouz[0]*By[6]-0.3535533905932737*By[4]*rhouz[5]+0.3535533905932737*Bz[4]*rhouy[5]+0.3535533905932737*Ex[4]*rho[5]+0.3535533905932737*rho[4]*Ex[5]+0.3535533905932737*rhouy[4]*Bz[5]-0.3535533905932737*rhouz[4]*By[5]-0.3535533905932737*By[2]*rhouz[3]+0.3535533905932737*Bz[2]*rhouy[3]+0.3535533905932737*Ex[2]*rho[3]+0.3535533905932737*rho[2]*Ex[3]+0.3535533905932737*rhouy[2]*Bz[3]-0.3535533905932737*rhouz[2]*By[3]; 
  outrhoux[7] += (-0.3535533905932737*By[0]*rhouz[7])+0.3535533905932737*Bz[0]*rhouy[7]+0.3535533905932737*Ex[0]*rho[7]+0.3535533905932737*rho[0]*Ex[7]+0.3535533905932737*rhouy[0]*Bz[7]-0.3535533905932737*rhouz[0]*By[7]-0.3535533905932737*By[1]*rhouz[6]+0.3535533905932737*Bz[1]*rhouy[6]+0.3535533905932737*Ex[1]*rho[6]+0.3535533905932737*rho[1]*Ex[6]+0.3535533905932737*rhouy[1]*Bz[6]-0.3535533905932737*rhouz[1]*By[6]-0.3535533905932737*By[2]*rhouz[5]+0.3535533905932737*Bz[2]*rhouy[5]+0.3535533905932737*Ex[2]*rho[5]+0.3535533905932737*rho[2]*Ex[5]+0.3535533905932737*rhouy[2]*Bz[5]-0.3535533905932737*rhouz[2]*By[5]-0.3535533905932737*By[3]*rhouz[4]+0.3535533905932737*Bz[3]*rhouy[4]+0.3535533905932737*Ex[3]*rho[4]+0.3535533905932737*rho[3]*Ex[4]+0.3535533905932737*rhouy[3]*Bz[4]-0.3535533905932737*rhouz[3]*By[4]; 

  outrhouy[0] += 0.3535533905932737*Bx[7]*rhouz[7]-0.3535533905932737*Bz[7]*rhoux[7]+0.3535533905932737*Ey[7]*rho[7]+0.3535533905932737*Bx[6]*rhouz[6]-0.3535533905932737*Bz[6]*rhoux[6]+0.3535533905932737*Ey[6]*rho[6]+0.3535533905932737*Bx[5]*rhouz[5]-0.3535533905932737*Bz[5]*rhoux[5]+0.3535533905932737*Ey[5]*rho[5]+0.3535533905932737*Bx[4]*rhouz[4]-0.3535533905932737*Bz[4]*rhoux[4]+0.3535533905932737*Ey[4]*rho[4]+0.3535533905932737*Bx[3]*rhouz[3]-0.3535533905932737*Bz[3]*rhoux[3]+0.3535533905932737*Ey[3]*rho[3]+0.3535533905932737*Bx[2]*rhouz[2]-0.3535533905932737*Bz[2]*rhoux[2]+0.3535533905932737*Ey[2]*rho[2]+0.3535533905932737*Bx[1]*rhouz[1]-0.3535533905932737*Bz[1]*rhoux[1]+0.3535533905932737*Ey[1]*rho[1]+0.3535533905932737*Bx[0]*rhouz[0]-0.3535533905932737*Bz[0]*rhoux[0]+0.3535533905932737*Ey[0]*rho[0]; 
  outrhouy[1] += 0.3535533905932737*Bx[6]*rhouz[7]-0.3535533905932737*Bz[6]*rhoux[7]+0.3535533905932737*Ey[6]*rho[7]+0.3535533905932737*rho[6]*Ey[7]-0.3535533905932737*rhoux[6]*Bz[7]+0.3535533905932737*rhouz[6]*Bx[7]+0.3535533905932737*Bx[3]*rhouz[5]-0.3535533905932737*Bz[3]*rhoux[5]+0.3535533905932737*Ey[3]*rho[5]+0.3535533905932737*rho[3]*Ey[5]-0.3535533905932737*rhoux[3]*Bz[5]+0.3535533905932737*rhouz[3]*Bx[5]+0.3535533905932737*Bx[2]*rhouz[4]-0.3535533905932737*Bz[2]*rhoux[4]+0.3535533905932737*Ey[2]*rho[4]+0.3535533905932737*rho[2]*Ey[4]-0.3535533905932737*rhoux[2]*Bz[4]+0.3535533905932737*rhouz[2]*Bx[4]+0.3535533905932737*Bx[0]*rhouz[1]-0.3535533905932737*Bz[0]*rhoux[1]+0.3535533905932737*Ey[0]*rho[1]+0.3535533905932737*rho[0]*Ey[1]-0.3535533905932737*rhoux[0]*Bz[1]+0.3535533905932737*rhouz[0]*Bx[1]; 
  outrhouy[2] += 0.3535533905932737*Bx[5]*rhouz[7]-0.3535533905932737*Bz[5]*rhoux[7]+0.3535533905932737*Ey[5]*rho[7]+0.3535533905932737*rho[5]*Ey[7]-0.3535533905932737*rhoux[5]*Bz[7]+0.3535533905932737*rhouz[5]*Bx[7]+0.3535533905932737*Bx[3]*rhouz[6]-0.3535533905932737*Bz[3]*rhoux[6]+0.3535533905932737*Ey[3]*rho[6]+0.3535533905932737*rho[3]*Ey[6]-0.3535533905932737*rhoux[3]*Bz[6]+0.3535533905932737*rhouz[3]*Bx[6]+0.3535533905932737*Bx[1]*rhouz[4]-0.3535533905932737*Bz[1]*rhoux[4]+0.3535533905932737*Ey[1]*rho[4]+0.3535533905932737*rho[1]*Ey[4]-0.3535533905932737*rhoux[1]*Bz[4]+0.3535533905932737*rhouz[1]*Bx[4]+0.3535533905932737*Bx[0]*rhouz[2]-0.3535533905932737*Bz[0]*rhoux[2]+0.3535533905932737*Ey[0]*rho[2]+0.3535533905932737*rho[0]*Ey[2]-0.3535533905932737*rhoux[0]*Bz[2]+0.3535533905932737*rhouz[0]*Bx[2]; 
  outrhouy[3] += 0.3535533905932737*Bx[4]*rhouz[7]-0.3535533905932737*Bz[4]*rhoux[7]+0.3535533905932737*Ey[4]*rho[7]+0.3535533905932737*rho[4]*Ey[7]-0.3535533905932737*rhoux[4]*Bz[7]+0.3535533905932737*rhouz[4]*Bx[7]+0.3535533905932737*Bx[2]*rhouz[6]-0.3535533905932737*Bz[2]*rhoux[6]+0.3535533905932737*Ey[2]*rho[6]+0.3535533905932737*rho[2]*Ey[6]-0.3535533905932737*rhoux[2]*Bz[6]+0.3535533905932737*rhouz[2]*Bx[6]+0.3535533905932737*Bx[1]*rhouz[5]-0.3535533905932737*Bz[1]*rhoux[5]+0.3535533905932737*Ey[1]*rho[5]+0.3535533905932737*rho[1]*Ey[5]-0.3535533905932737*rhoux[1]*Bz[5]+0.3535533905932737*rhouz[1]*Bx[5]+0.3535533905932737*Bx[0]*rhouz[3]-0.3535533905932737*Bz[0]*rhoux[3]+0.3535533905932737*Ey[0]*rho[3]+0.3535533905932737*rho[0]*Ey[3]-0.3535533905932737*rhoux[0]*Bz[3]+0.3535533905932737*rhouz[0]*Bx[3]; 
  outrhouy[4] += 0.3535533905932737*Bx[3]*rhouz[7]-0.3535533905932737*Bz[3]*rhoux[7]+0.3535533905932737*Ey[3]*rho[7]+0.3535533905932737*rho[3]*Ey[7]-0.3535533905932737*rhoux[3]*Bz[7]+0.3535533905932737*rhouz[3]*Bx[7]+0.3535533905932737*Bx[5]*rhouz[6]-0.3535533905932737*Bz[5]*rhoux[6]+0.3535533905932737*Ey[5]*rho[6]+0.3535533905932737*rho[5]*Ey[6]-0.3535533905932737*rhoux[5]*Bz[6]+0.3535533905932737*rhouz[5]*Bx[6]+0.3535533905932737*Bx[0]*rhouz[4]-0.3535533905932737*Bz[0]*rhoux[4]+0.3535533905932737*Ey[0]*rho[4]+0.3535533905932737*rho[0]*Ey[4]-0.3535533905932737*rhoux[0]*Bz[4]+0.3535533905932737*rhouz[0]*Bx[4]+0.3535533905932737*Bx[1]*rhouz[2]-0.3535533905932737*Bz[1]*rhoux[2]+0.3535533905932737*Ey[1]*rho[2]+0.3535533905932737*rho[1]*Ey[2]-0.3535533905932737*rhoux[1]*Bz[2]+0.3535533905932737*rhouz[1]*Bx[2]; 
  outrhouy[5] += 0.3535533905932737*Bx[2]*rhouz[7]-0.3535533905932737*Bz[2]*rhoux[7]+0.3535533905932737*Ey[2]*rho[7]+0.3535533905932737*rho[2]*Ey[7]-0.3535533905932737*rhoux[2]*Bz[7]+0.3535533905932737*rhouz[2]*Bx[7]+0.3535533905932737*Bx[4]*rhouz[6]-0.3535533905932737*Bz[4]*rhoux[6]+0.3535533905932737*Ey[4]*rho[6]+0.3535533905932737*rho[4]*Ey[6]-0.3535533905932737*rhoux[4]*Bz[6]+0.3535533905932737*rhouz[4]*Bx[6]+0.3535533905932737*Bx[0]*rhouz[5]-0.3535533905932737*Bz[0]*rhoux[5]+0.3535533905932737*Ey[0]*rho[5]+0.3535533905932737*rho[0]*Ey[5]-0.3535533905932737*rhoux[0]*Bz[5]+0.3535533905932737*rhouz[0]*Bx[5]+0.3535533905932737*Bx[1]*rhouz[3]-0.3535533905932737*Bz[1]*rhoux[3]+0.3535533905932737*Ey[1]*rho[3]+0.3535533905932737*rho[1]*Ey[3]-0.3535533905932737*rhoux[1]*Bz[3]+0.3535533905932737*rhouz[1]*Bx[3]; 
  outrhouy[6] += 0.3535533905932737*Bx[1]*rhouz[7]-0.3535533905932737*Bz[1]*rhoux[7]+0.3535533905932737*Ey[1]*rho[7]+0.3535533905932737*rho[1]*Ey[7]-0.3535533905932737*rhoux[1]*Bz[7]+0.3535533905932737*rhouz[1]*Bx[7]+0.3535533905932737*Bx[0]*rhouz[6]-0.3535533905932737*Bz[0]*rhoux[6]+0.3535533905932737*Ey[0]*rho[6]+0.3535533905932737*rho[0]*Ey[6]-0.3535533905932737*rhoux[0]*Bz[6]+0.3535533905932737*rhouz[0]*Bx[6]+0.3535533905932737*Bx[4]*rhouz[5]-0.3535533905932737*Bz[4]*rhoux[5]+0.3535533905932737*Ey[4]*rho[5]+0.3535533905932737*rho[4]*Ey[5]-0.3535533905932737*rhoux[4]*Bz[5]+0.3535533905932737*rhouz[4]*Bx[5]+0.3535533905932737*Bx[2]*rhouz[3]-0.3535533905932737*Bz[2]*rhoux[3]+0.3535533905932737*Ey[2]*rho[3]+0.3535533905932737*rho[2]*Ey[3]-0.3535533905932737*rhoux[2]*Bz[3]+0.3535533905932737*rhouz[2]*Bx[3]; 
  outrhouy[7] += 0.3535533905932737*Bx[0]*rhouz[7]-0.3535533905932737*Bz[0]*rhoux[7]+0.3535533905932737*Ey[0]*rho[7]+0.3535533905932737*rho[0]*Ey[7]-0.3535533905932737*rhoux[0]*Bz[7]+0.3535533905932737*rhouz[0]*Bx[7]+0.3535533905932737*Bx[1]*rhouz[6]-0.3535533905932737*Bz[1]*rhoux[6]+0.3535533905932737*Ey[1]*rho[6]+0.3535533905932737*rho[1]*Ey[6]-0.3535533905932737*rhoux[1]*Bz[6]+0.3535533905932737*rhouz[1]*Bx[6]+0.3535533905932737*Bx[2]*rhouz[5]-0.3535533905932737*Bz[2]*rhoux[5]+0.3535533905932737*Ey[2]*rho[5]+0.3535533905932737*rho[2]*Ey[5]-0.3535533905932737*rhoux[2]*Bz[5]+0.3535533905932737*rhouz[2]*Bx[5]+0.3535533905932737*Bx[3]*rhouz[4]-0.3535533905932737*Bz[3]*rhoux[4]+0.3535533905932737*Ey[3]*rho[4]+0.3535533905932737*rho[3]*Ey[4]-0.3535533905932737*rhoux[3]*Bz[4]+0.3535533905932737*rhouz[3]*Bx[4]; 

  outrhouz[0] += (-0.3535533905932737*Bx[7]*rhouy[7])+0.3535533905932737*By[7]*rhoux[7]+0.3535533905932737*Ez[7]*rho[7]-0.3535533905932737*Bx[6]*rhouy[6]+0.3535533905932737*By[6]*rhoux[6]+0.3535533905932737*Ez[6]*rho[6]-0.3535533905932737*Bx[5]*rhouy[5]+0.3535533905932737*By[5]*rhoux[5]+0.3535533905932737*Ez[5]*rho[5]-0.3535533905932737*Bx[4]*rhouy[4]+0.3535533905932737*By[4]*rhoux[4]+0.3535533905932737*Ez[4]*rho[4]-0.3535533905932737*Bx[3]*rhouy[3]+0.3535533905932737*By[3]*rhoux[3]+0.3535533905932737*Ez[3]*rho[3]-0.3535533905932737*Bx[2]*rhouy[2]+0.3535533905932737*By[2]*rhoux[2]+0.3535533905932737*Ez[2]*rho[2]-0.3535533905932737*Bx[1]*rhouy[1]+0.3535533905932737*By[1]*rhoux[1]+0.3535533905932737*Ez[1]*rho[1]-0.3535533905932737*Bx[0]*rhouy[0]+0.3535533905932737*By[0]*rhoux[0]+0.3535533905932737*Ez[0]*rho[0]; 
  outrhouz[1] += (-0.3535533905932737*Bx[6]*rhouy[7])+0.3535533905932737*By[6]*rhoux[7]+0.3535533905932737*Ez[6]*rho[7]+0.3535533905932737*rho[6]*Ez[7]+0.3535533905932737*rhoux[6]*By[7]-0.3535533905932737*rhouy[6]*Bx[7]-0.3535533905932737*Bx[3]*rhouy[5]+0.3535533905932737*By[3]*rhoux[5]+0.3535533905932737*Ez[3]*rho[5]+0.3535533905932737*rho[3]*Ez[5]+0.3535533905932737*rhoux[3]*By[5]-0.3535533905932737*rhouy[3]*Bx[5]-0.3535533905932737*Bx[2]*rhouy[4]+0.3535533905932737*By[2]*rhoux[4]+0.3535533905932737*Ez[2]*rho[4]+0.3535533905932737*rho[2]*Ez[4]+0.3535533905932737*rhoux[2]*By[4]-0.3535533905932737*rhouy[2]*Bx[4]-0.3535533905932737*Bx[0]*rhouy[1]+0.3535533905932737*By[0]*rhoux[1]+0.3535533905932737*Ez[0]*rho[1]+0.3535533905932737*rho[0]*Ez[1]+0.3535533905932737*rhoux[0]*By[1]-0.3535533905932737*rhouy[0]*Bx[1]; 
  outrhouz[2] += (-0.3535533905932737*Bx[5]*rhouy[7])+0.3535533905932737*By[5]*rhoux[7]+0.3535533905932737*Ez[5]*rho[7]+0.3535533905932737*rho[5]*Ez[7]+0.3535533905932737*rhoux[5]*By[7]-0.3535533905932737*rhouy[5]*Bx[7]-0.3535533905932737*Bx[3]*rhouy[6]+0.3535533905932737*By[3]*rhoux[6]+0.3535533905932737*Ez[3]*rho[6]+0.3535533905932737*rho[3]*Ez[6]+0.3535533905932737*rhoux[3]*By[6]-0.3535533905932737*rhouy[3]*Bx[6]-0.3535533905932737*Bx[1]*rhouy[4]+0.3535533905932737*By[1]*rhoux[4]+0.3535533905932737*Ez[1]*rho[4]+0.3535533905932737*rho[1]*Ez[4]+0.3535533905932737*rhoux[1]*By[4]-0.3535533905932737*rhouy[1]*Bx[4]-0.3535533905932737*Bx[0]*rhouy[2]+0.3535533905932737*By[0]*rhoux[2]+0.3535533905932737*Ez[0]*rho[2]+0.3535533905932737*rho[0]*Ez[2]+0.3535533905932737*rhoux[0]*By[2]-0.3535533905932737*rhouy[0]*Bx[2]; 
  outrhouz[3] += (-0.3535533905932737*Bx[4]*rhouy[7])+0.3535533905932737*By[4]*rhoux[7]+0.3535533905932737*Ez[4]*rho[7]+0.3535533905932737*rho[4]*Ez[7]+0.3535533905932737*rhoux[4]*By[7]-0.3535533905932737*rhouy[4]*Bx[7]-0.3535533905932737*Bx[2]*rhouy[6]+0.3535533905932737*By[2]*rhoux[6]+0.3535533905932737*Ez[2]*rho[6]+0.3535533905932737*rho[2]*Ez[6]+0.3535533905932737*rhoux[2]*By[6]-0.3535533905932737*rhouy[2]*Bx[6]-0.3535533905932737*Bx[1]*rhouy[5]+0.3535533905932737*By[1]*rhoux[5]+0.3535533905932737*Ez[1]*rho[5]+0.3535533905932737*rho[1]*Ez[5]+0.3535533905932737*rhoux[1]*By[5]-0.3535533905932737*rhouy[1]*Bx[5]-0.3535533905932737*Bx[0]*rhouy[3]+0.3535533905932737*By[0]*rhoux[3]+0.3535533905932737*Ez[0]*rho[3]+0.3535533905932737*rho[0]*Ez[3]+0.3535533905932737*rhoux[0]*By[3]-0.3535533905932737*rhouy[0]*Bx[3]; 
  outrhouz[4] += (-0.3535533905932737*Bx[3]*rhouy[7])+0.3535533905932737*By[3]*rhoux[7]+0.3535533905932737*Ez[3]*rho[7]+0.3535533905932737*rho[3]*Ez[7]+0.3535533905932737*rhoux[3]*By[7]-0.3535533905932737*rhouy[3]*Bx[7]-0.3535533905932737*Bx[5]*rhouy[6]+0.3535533905932737*By[5]*rhoux[6]+0.3535533905932737*Ez[5]*rho[6]+0.3535533905932737*rho[5]*Ez[6]+0.3535533905932737*rhoux[5]*By[6]-0.3535533905932737*rhouy[5]*Bx[6]-0.3535533905932737*Bx[0]*rhouy[4]+0.3535533905932737*By[0]*rhoux[4]+0.3535533905932737*Ez[0]*rho[4]+0.3535533905932737*rho[0]*Ez[4]+0.3535533905932737*rhoux[0]*By[4]-0.3535533905932737*rhouy[0]*Bx[4]-0.3535533905932737*Bx[1]*rhouy[2]+0.3535533905932737*By[1]*rhoux[2]+0.3535533905932737*Ez[1]*rho[2]+0.3535533905932737*rho[1]*Ez[2]+0.3535533905932737*rhoux[1]*By[2]-0.3535533905932737*rhouy[1]*Bx[2]; 
  outrhouz[5] += (-0.3535533905932737*Bx[2]*rhouy[7])+0.3535533905932737*By[2]*rhoux[7]+0.3535533905932737*Ez[2]*rho[7]+0.3535533905932737*rho[2]*Ez[7]+0.3535533905932737*rhoux[2]*By[7]-0.3535533905932737*rhouy[2]*Bx[7]-0.3535533905932737*Bx[4]*rhouy[6]+0.3535533905932737*By[4]*rhoux[6]+0.3535533905932737*Ez[4]*rho[6]+0.3535533905932737*rho[4]*Ez[6]+0.3535533905932737*rhoux[4]*By[6]-0.3535533905932737*rhouy[4]*Bx[6]-0.3535533905932737*Bx[0]*rhouy[5]+0.3535533905932737*By[0]*rhoux[5]+0.3535533905932737*Ez[0]*rho[5]+0.3535533905932737*rho[0]*Ez[5]+0.3535533905932737*rhoux[0]*By[5]-0.3535533905932737*rhouy[0]*Bx[5]-0.3535533905932737*Bx[1]*rhouy[3]+0.3535533905932737*By[1]*rhoux[3]+0.3535533905932737*Ez[1]*rho[3]+0.3535533905932737*rho[1]*Ez[3]+0.3535533905932737*rhoux[1]*By[3]-0.3535533905932737*rhouy[1]*Bx[3]; 
  outrhouz[6] += (-0.3535533905932737*Bx[1]*rhouy[7])+0.3535533905932737*By[1]*rhoux[7]+0.3535533905932737*Ez[1]*rho[7]+0.3535533905932737*rho[1]*Ez[7]+0.3535533905932737*rhoux[1]*By[7]-0.3535533905932737*rhouy[1]*Bx[7]-0.3535533905932737*Bx[0]*rhouy[6]+0.3535533905932737*By[0]*rhoux[6]+0.3535533905932737*Ez[0]*rho[6]+0.3535533905932737*rho[0]*Ez[6]+0.3535533905932737*rhoux[0]*By[6]-0.3535533905932737*rhouy[0]*Bx[6]-0.3535533905932737*Bx[4]*rhouy[5]+0.3535533905932737*By[4]*rhoux[5]+0.3535533905932737*Ez[4]*rho[5]+0.3535533905932737*rho[4]*Ez[5]+0.3535533905932737*rhoux[4]*By[5]-0.3535533905932737*rhouy[4]*Bx[5]-0.3535533905932737*Bx[2]*rhouy[3]+0.3535533905932737*By[2]*rhoux[3]+0.3535533905932737*Ez[2]*rho[3]+0.3535533905932737*rho[2]*Ez[3]+0.3535533905932737*rhoux[2]*By[3]-0.3535533905932737*rhouy[2]*Bx[3]; 
  outrhouz[7] += (-0.3535533905932737*Bx[0]*rhouy[7])+0.3535533905932737*By[0]*rhoux[7]+0.3535533905932737*Ez[0]*rho[7]+0.3535533905932737*rho[0]*Ez[7]+0.3535533905932737*rhoux[0]*By[7]-0.3535533905932737*rhouy[0]*Bx[7]-0.3535533905932737*Bx[1]*rhouy[6]+0.3535533905932737*By[1]*rhoux[6]+0.3535533905932737*Ez[1]*rho[6]+0.3535533905932737*rho[1]*Ez[6]+0.3535533905932737*rhoux[1]*By[6]-0.3535533905932737*rhouy[1]*Bx[6]-0.3535533905932737*Bx[2]*rhouy[5]+0.3535533905932737*By[2]*rhoux[5]+0.3535533905932737*Ez[2]*rho[5]+0.3535533905932737*rho[2]*Ez[5]+0.3535533905932737*rhoux[2]*By[5]-0.3535533905932737*rhouy[2]*Bx[5]-0.3535533905932737*Bx[3]*rhouy[4]+0.3535533905932737*By[3]*rhoux[4]+0.3535533905932737*Ez[3]*rho[4]+0.3535533905932737*rho[3]*Ez[4]+0.3535533905932737*rhoux[3]*By[4]-0.3535533905932737*rhouy[3]*Bx[4]; 

} 
