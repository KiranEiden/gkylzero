#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_source_2x_ser_p2(const double *qmem, const double *nu, const double *nu_vth_sq, const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *rhou_perp_i, const double *p_perp, double* GKYL_RESTRICT out) 
{ 
  // qmem:             q/m*EM fields.
  // nu:               Collisionality.
  // nu_vth_sq:        nu*vth^2 = nu*T/m.
  // vlasov_pkpm_moms: [rho, p_parallel, q_parallel], Moments computed from kinetic equation in pkpm model.
  // euler_pkpm:       Input fluid variables.
  // rhou_perp_i:      Input perpendicular momentum density, rhou - (rhou . b)b = [rhou_perp_x, rhou_perp_y, rhou_perp_z].
  // p_perp:           Input perpendicular pressure (E_perp = 1/2 rho u_perp^2 + p_perp).
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

  const double *rhou_perp_x = &rhou_perp_i[0]; 
  const double *rhou_perp_y = &rhou_perp_i[8]; 
  const double *rhou_perp_z = &rhou_perp_i[16]; 

  double *outrhoux = &out[0]; 
  double *outrhouy = &out[8]; 
  double *outrhouz = &out[16]; 
  double *outE_perp = &out[24]; 

  outrhoux[0] += (-0.5*By[7]*rhouz[7])+0.5*Bz[7]*rhouy[7]+0.5*Ex[7]*rho[7]-0.5*By[6]*rhouz[6]+0.5*Bz[6]*rhouy[6]+0.5*Ex[6]*rho[6]-0.5*By[5]*rhouz[5]+0.5*Bz[5]*rhouy[5]+0.5*Ex[5]*rho[5]-0.5*By[4]*rhouz[4]+0.5*Bz[4]*rhouy[4]+0.5*Ex[4]*rho[4]-0.5*By[3]*rhouz[3]+0.5*Bz[3]*rhouy[3]+0.5*Ex[3]*rho[3]-0.5*By[2]*rhouz[2]+0.5*Bz[2]*rhouy[2]+0.5*Ex[2]*rho[2]-0.5*By[1]*rhouz[1]+0.5*Bz[1]*rhouy[1]+0.5*Ex[1]*rho[1]-0.5*By[0]*rhouz[0]+0.5*Bz[0]*rhouy[0]+0.5*Ex[0]*rho[0]; 
  outrhoux[1] += (-0.5000000000000001*By[5]*rhouz[7])+0.5000000000000001*Bz[5]*rhouy[7]+0.5000000000000001*Ex[5]*rho[7]+0.5000000000000001*rho[5]*Ex[7]+0.5000000000000001*rhouy[5]*Bz[7]-0.5000000000000001*rhouz[5]*By[7]-0.447213595499958*By[3]*rhouz[6]+0.447213595499958*Bz[3]*rhouy[6]+0.447213595499958*Ex[3]*rho[6]+0.447213595499958*rho[3]*Ex[6]+0.447213595499958*rhouy[3]*Bz[6]-0.447213595499958*rhouz[3]*By[6]-0.4472135954999579*By[1]*rhouz[4]+0.4472135954999579*Bz[1]*rhouy[4]+0.4472135954999579*Ex[1]*rho[4]+0.4472135954999579*rho[1]*Ex[4]+0.4472135954999579*rhouy[1]*Bz[4]-0.4472135954999579*rhouz[1]*By[4]-0.5*By[2]*rhouz[3]+0.5*Bz[2]*rhouy[3]+0.5*Ex[2]*rho[3]+0.5*rho[2]*Ex[3]+0.5*rhouy[2]*Bz[3]-0.5*rhouz[2]*By[3]-0.5*By[0]*rhouz[1]+0.5*Bz[0]*rhouy[1]+0.5*Ex[0]*rho[1]+0.5*rho[0]*Ex[1]+0.5*rhouy[0]*Bz[1]-0.5*rhouz[0]*By[1]; 
  outrhoux[2] += (-0.447213595499958*By[3]*rhouz[7])+0.447213595499958*Bz[3]*rhouy[7]+0.447213595499958*Ex[3]*rho[7]+0.447213595499958*rho[3]*Ex[7]+0.447213595499958*rhouy[3]*Bz[7]-0.447213595499958*rhouz[3]*By[7]-0.5000000000000001*By[4]*rhouz[6]+0.5000000000000001*Bz[4]*rhouy[6]+0.5000000000000001*Ex[4]*rho[6]+0.5000000000000001*rho[4]*Ex[6]+0.5000000000000001*rhouy[4]*Bz[6]-0.5000000000000001*rhouz[4]*By[6]-0.4472135954999579*By[2]*rhouz[5]+0.4472135954999579*Bz[2]*rhouy[5]+0.4472135954999579*Ex[2]*rho[5]+0.4472135954999579*rho[2]*Ex[5]+0.4472135954999579*rhouy[2]*Bz[5]-0.4472135954999579*rhouz[2]*By[5]-0.5*By[1]*rhouz[3]+0.5*Bz[1]*rhouy[3]+0.5*Ex[1]*rho[3]+0.5*rho[1]*Ex[3]+0.5*rhouy[1]*Bz[3]-0.5*rhouz[1]*By[3]-0.5*By[0]*rhouz[2]+0.5*Bz[0]*rhouy[2]+0.5*Ex[0]*rho[2]+0.5*rho[0]*Ex[2]+0.5*rhouy[0]*Bz[2]-0.5*rhouz[0]*By[2]; 
  outrhoux[3] += (-0.4*By[6]*rhouz[7])-0.447213595499958*By[2]*rhouz[7]+0.4*Bz[6]*rhouy[7]+0.447213595499958*Bz[2]*rhouy[7]+0.4*Ex[6]*rho[7]+0.447213595499958*Ex[2]*rho[7]+0.4*rho[6]*Ex[7]+0.447213595499958*rho[2]*Ex[7]+0.4*rhouy[6]*Bz[7]+0.447213595499958*rhouy[2]*Bz[7]-0.4*rhouz[6]*By[7]-0.447213595499958*rhouz[2]*By[7]-0.447213595499958*By[1]*rhouz[6]+0.447213595499958*Bz[1]*rhouy[6]+0.447213595499958*Ex[1]*rho[6]+0.447213595499958*rho[1]*Ex[6]+0.447213595499958*rhouy[1]*Bz[6]-0.447213595499958*rhouz[1]*By[6]-0.4472135954999579*By[3]*rhouz[5]+0.4472135954999579*Bz[3]*rhouy[5]+0.4472135954999579*Ex[3]*rho[5]+0.4472135954999579*rho[3]*Ex[5]+0.4472135954999579*rhouy[3]*Bz[5]-0.4472135954999579*rhouz[3]*By[5]-0.4472135954999579*By[3]*rhouz[4]+0.4472135954999579*Bz[3]*rhouy[4]+0.4472135954999579*Ex[3]*rho[4]+0.4472135954999579*rho[3]*Ex[4]+0.4472135954999579*rhouy[3]*Bz[4]-0.4472135954999579*rhouz[3]*By[4]-0.5*By[0]*rhouz[3]+0.5*Bz[0]*rhouy[3]+0.5*Ex[0]*rho[3]+0.5*rho[0]*Ex[3]+0.5*rhouy[0]*Bz[3]-0.5*rhouz[0]*By[3]-0.5*By[1]*rhouz[2]+0.5*Bz[1]*rhouy[2]+0.5*Ex[1]*rho[2]+0.5*rho[1]*Ex[2]+0.5*rhouy[1]*Bz[2]-0.5*rhouz[1]*By[2]; 
  outrhoux[4] += (-0.4472135954999579*By[7]*rhouz[7])+0.4472135954999579*Bz[7]*rhouy[7]+0.4472135954999579*Ex[7]*rho[7]-0.31943828249997*By[6]*rhouz[6]-0.5000000000000001*By[2]*rhouz[6]+0.31943828249997*Bz[6]*rhouy[6]+0.5000000000000001*Bz[2]*rhouy[6]+0.31943828249997*Ex[6]*rho[6]+0.5000000000000001*Ex[2]*rho[6]+0.5000000000000001*rho[2]*Ex[6]+0.5000000000000001*rhouy[2]*Bz[6]-0.5000000000000001*rhouz[2]*By[6]-0.31943828249997*By[4]*rhouz[4]-0.5*By[0]*rhouz[4]+0.31943828249997*Bz[4]*rhouy[4]+0.5*Bz[0]*rhouy[4]+0.31943828249997*Ex[4]*rho[4]+0.5*Ex[0]*rho[4]+0.5*rho[0]*Ex[4]+0.5*rhouy[0]*Bz[4]-0.5*rhouz[0]*By[4]-0.4472135954999579*By[3]*rhouz[3]+0.4472135954999579*Bz[3]*rhouy[3]+0.4472135954999579*Ex[3]*rho[3]-0.4472135954999579*By[1]*rhouz[1]+0.4472135954999579*Bz[1]*rhouy[1]+0.4472135954999579*Ex[1]*rho[1]; 
  outrhoux[5] += (-0.31943828249997*By[7]*rhouz[7])-0.5000000000000001*By[1]*rhouz[7]+0.31943828249997*Bz[7]*rhouy[7]+0.5000000000000001*Bz[1]*rhouy[7]+0.31943828249997*Ex[7]*rho[7]+0.5000000000000001*Ex[1]*rho[7]+0.5000000000000001*rho[1]*Ex[7]+0.5000000000000001*rhouy[1]*Bz[7]-0.5000000000000001*rhouz[1]*By[7]-0.4472135954999579*By[6]*rhouz[6]+0.4472135954999579*Bz[6]*rhouy[6]+0.4472135954999579*Ex[6]*rho[6]-0.31943828249997*By[5]*rhouz[5]-0.5*By[0]*rhouz[5]+0.31943828249997*Bz[5]*rhouy[5]+0.5*Bz[0]*rhouy[5]+0.31943828249997*Ex[5]*rho[5]+0.5*Ex[0]*rho[5]+0.5*rho[0]*Ex[5]+0.5*rhouy[0]*Bz[5]-0.5*rhouz[0]*By[5]-0.4472135954999579*By[3]*rhouz[3]+0.4472135954999579*Bz[3]*rhouy[3]+0.4472135954999579*Ex[3]*rho[3]-0.4472135954999579*By[2]*rhouz[2]+0.4472135954999579*Bz[2]*rhouy[2]+0.4472135954999579*Ex[2]*rho[2]; 
  outrhoux[6] += (-0.4*By[3]*rhouz[7])+0.4*Bz[3]*rhouy[7]+0.4*Ex[3]*rho[7]+0.4*rho[3]*Ex[7]+0.4*rhouy[3]*Bz[7]-0.4*rhouz[3]*By[7]-0.4472135954999579*By[5]*rhouz[6]-0.31943828249997*By[4]*rhouz[6]-0.5*By[0]*rhouz[6]+0.4472135954999579*Bz[5]*rhouy[6]+0.31943828249997*Bz[4]*rhouy[6]+0.5*Bz[0]*rhouy[6]+0.4472135954999579*Ex[5]*rho[6]+0.31943828249997*Ex[4]*rho[6]+0.5*Ex[0]*rho[6]+0.4472135954999579*rho[5]*Ex[6]+0.31943828249997*rho[4]*Ex[6]+0.5*rho[0]*Ex[6]+0.4472135954999579*rhouy[5]*Bz[6]+0.31943828249997*rhouy[4]*Bz[6]+0.5*rhouy[0]*Bz[6]-0.4472135954999579*rhouz[5]*By[6]-0.31943828249997*rhouz[4]*By[6]-0.5*rhouz[0]*By[6]-0.5000000000000001*By[2]*rhouz[4]+0.5000000000000001*Bz[2]*rhouy[4]+0.5000000000000001*Ex[2]*rho[4]+0.5000000000000001*rho[2]*Ex[4]+0.5000000000000001*rhouy[2]*Bz[4]-0.5000000000000001*rhouz[2]*By[4]-0.447213595499958*By[1]*rhouz[3]+0.447213595499958*Bz[1]*rhouy[3]+0.447213595499958*Ex[1]*rho[3]+0.447213595499958*rho[1]*Ex[3]+0.447213595499958*rhouy[1]*Bz[3]-0.447213595499958*rhouz[1]*By[3]; 
  outrhoux[7] += (-0.31943828249997*By[5]*rhouz[7])-0.4472135954999579*By[4]*rhouz[7]-0.5*By[0]*rhouz[7]+0.31943828249997*Bz[5]*rhouy[7]+0.4472135954999579*Bz[4]*rhouy[7]+0.5*Bz[0]*rhouy[7]+0.31943828249997*Ex[5]*rho[7]+0.4472135954999579*Ex[4]*rho[7]+0.5*Ex[0]*rho[7]+0.31943828249997*rho[5]*Ex[7]+0.4472135954999579*rho[4]*Ex[7]+0.5*rho[0]*Ex[7]+0.31943828249997*rhouy[5]*Bz[7]+0.4472135954999579*rhouy[4]*Bz[7]+0.5*rhouy[0]*Bz[7]-0.31943828249997*rhouz[5]*By[7]-0.4472135954999579*rhouz[4]*By[7]-0.5*rhouz[0]*By[7]-0.4*By[3]*rhouz[6]+0.4*Bz[3]*rhouy[6]+0.4*Ex[3]*rho[6]+0.4*rho[3]*Ex[6]+0.4*rhouy[3]*Bz[6]-0.4*rhouz[3]*By[6]-0.5000000000000001*By[1]*rhouz[5]+0.5000000000000001*Bz[1]*rhouy[5]+0.5000000000000001*Ex[1]*rho[5]+0.5000000000000001*rho[1]*Ex[5]+0.5000000000000001*rhouy[1]*Bz[5]-0.5000000000000001*rhouz[1]*By[5]-0.447213595499958*By[2]*rhouz[3]+0.447213595499958*Bz[2]*rhouy[3]+0.447213595499958*Ex[2]*rho[3]+0.447213595499958*rho[2]*Ex[3]+0.447213595499958*rhouy[2]*Bz[3]-0.447213595499958*rhouz[2]*By[3]; 

  outrhouy[0] += 0.5*Bx[7]*rhouz[7]-0.5*Bz[7]*rhoux[7]+0.5*Ey[7]*rho[7]+0.5*Bx[6]*rhouz[6]-0.5*Bz[6]*rhoux[6]+0.5*Ey[6]*rho[6]+0.5*Bx[5]*rhouz[5]-0.5*Bz[5]*rhoux[5]+0.5*Ey[5]*rho[5]+0.5*Bx[4]*rhouz[4]-0.5*Bz[4]*rhoux[4]+0.5*Ey[4]*rho[4]+0.5*Bx[3]*rhouz[3]-0.5*Bz[3]*rhoux[3]+0.5*Ey[3]*rho[3]+0.5*Bx[2]*rhouz[2]-0.5*Bz[2]*rhoux[2]+0.5*Ey[2]*rho[2]+0.5*Bx[1]*rhouz[1]-0.5*Bz[1]*rhoux[1]+0.5*Ey[1]*rho[1]+0.5*Bx[0]*rhouz[0]-0.5*Bz[0]*rhoux[0]+0.5*Ey[0]*rho[0]; 
  outrhouy[1] += 0.5000000000000001*Bx[5]*rhouz[7]-0.5000000000000001*Bz[5]*rhoux[7]+0.5000000000000001*Ey[5]*rho[7]+0.5000000000000001*rho[5]*Ey[7]-0.5000000000000001*rhoux[5]*Bz[7]+0.5000000000000001*rhouz[5]*Bx[7]+0.447213595499958*Bx[3]*rhouz[6]-0.447213595499958*Bz[3]*rhoux[6]+0.447213595499958*Ey[3]*rho[6]+0.447213595499958*rho[3]*Ey[6]-0.447213595499958*rhoux[3]*Bz[6]+0.447213595499958*rhouz[3]*Bx[6]+0.4472135954999579*Bx[1]*rhouz[4]-0.4472135954999579*Bz[1]*rhoux[4]+0.4472135954999579*Ey[1]*rho[4]+0.4472135954999579*rho[1]*Ey[4]-0.4472135954999579*rhoux[1]*Bz[4]+0.4472135954999579*rhouz[1]*Bx[4]+0.5*Bx[2]*rhouz[3]-0.5*Bz[2]*rhoux[3]+0.5*Ey[2]*rho[3]+0.5*rho[2]*Ey[3]-0.5*rhoux[2]*Bz[3]+0.5*rhouz[2]*Bx[3]+0.5*Bx[0]*rhouz[1]-0.5*Bz[0]*rhoux[1]+0.5*Ey[0]*rho[1]+0.5*rho[0]*Ey[1]-0.5*rhoux[0]*Bz[1]+0.5*rhouz[0]*Bx[1]; 
  outrhouy[2] += 0.447213595499958*Bx[3]*rhouz[7]-0.447213595499958*Bz[3]*rhoux[7]+0.447213595499958*Ey[3]*rho[7]+0.447213595499958*rho[3]*Ey[7]-0.447213595499958*rhoux[3]*Bz[7]+0.447213595499958*rhouz[3]*Bx[7]+0.5000000000000001*Bx[4]*rhouz[6]-0.5000000000000001*Bz[4]*rhoux[6]+0.5000000000000001*Ey[4]*rho[6]+0.5000000000000001*rho[4]*Ey[6]-0.5000000000000001*rhoux[4]*Bz[6]+0.5000000000000001*rhouz[4]*Bx[6]+0.4472135954999579*Bx[2]*rhouz[5]-0.4472135954999579*Bz[2]*rhoux[5]+0.4472135954999579*Ey[2]*rho[5]+0.4472135954999579*rho[2]*Ey[5]-0.4472135954999579*rhoux[2]*Bz[5]+0.4472135954999579*rhouz[2]*Bx[5]+0.5*Bx[1]*rhouz[3]-0.5*Bz[1]*rhoux[3]+0.5*Ey[1]*rho[3]+0.5*rho[1]*Ey[3]-0.5*rhoux[1]*Bz[3]+0.5*rhouz[1]*Bx[3]+0.5*Bx[0]*rhouz[2]-0.5*Bz[0]*rhoux[2]+0.5*Ey[0]*rho[2]+0.5*rho[0]*Ey[2]-0.5*rhoux[0]*Bz[2]+0.5*rhouz[0]*Bx[2]; 
  outrhouy[3] += 0.4*Bx[6]*rhouz[7]+0.447213595499958*Bx[2]*rhouz[7]-0.4*Bz[6]*rhoux[7]-0.447213595499958*Bz[2]*rhoux[7]+0.4*Ey[6]*rho[7]+0.447213595499958*Ey[2]*rho[7]+0.4*rho[6]*Ey[7]+0.447213595499958*rho[2]*Ey[7]-0.4*rhoux[6]*Bz[7]-0.447213595499958*rhoux[2]*Bz[7]+0.4*rhouz[6]*Bx[7]+0.447213595499958*rhouz[2]*Bx[7]+0.447213595499958*Bx[1]*rhouz[6]-0.447213595499958*Bz[1]*rhoux[6]+0.447213595499958*Ey[1]*rho[6]+0.447213595499958*rho[1]*Ey[6]-0.447213595499958*rhoux[1]*Bz[6]+0.447213595499958*rhouz[1]*Bx[6]+0.4472135954999579*Bx[3]*rhouz[5]-0.4472135954999579*Bz[3]*rhoux[5]+0.4472135954999579*Ey[3]*rho[5]+0.4472135954999579*rho[3]*Ey[5]-0.4472135954999579*rhoux[3]*Bz[5]+0.4472135954999579*rhouz[3]*Bx[5]+0.4472135954999579*Bx[3]*rhouz[4]-0.4472135954999579*Bz[3]*rhoux[4]+0.4472135954999579*Ey[3]*rho[4]+0.4472135954999579*rho[3]*Ey[4]-0.4472135954999579*rhoux[3]*Bz[4]+0.4472135954999579*rhouz[3]*Bx[4]+0.5*Bx[0]*rhouz[3]-0.5*Bz[0]*rhoux[3]+0.5*Ey[0]*rho[3]+0.5*rho[0]*Ey[3]-0.5*rhoux[0]*Bz[3]+0.5*rhouz[0]*Bx[3]+0.5*Bx[1]*rhouz[2]-0.5*Bz[1]*rhoux[2]+0.5*Ey[1]*rho[2]+0.5*rho[1]*Ey[2]-0.5*rhoux[1]*Bz[2]+0.5*rhouz[1]*Bx[2]; 
  outrhouy[4] += 0.4472135954999579*Bx[7]*rhouz[7]-0.4472135954999579*Bz[7]*rhoux[7]+0.4472135954999579*Ey[7]*rho[7]+0.31943828249997*Bx[6]*rhouz[6]+0.5000000000000001*Bx[2]*rhouz[6]-0.31943828249997*Bz[6]*rhoux[6]-0.5000000000000001*Bz[2]*rhoux[6]+0.31943828249997*Ey[6]*rho[6]+0.5000000000000001*Ey[2]*rho[6]+0.5000000000000001*rho[2]*Ey[6]-0.5000000000000001*rhoux[2]*Bz[6]+0.5000000000000001*rhouz[2]*Bx[6]+0.31943828249997*Bx[4]*rhouz[4]+0.5*Bx[0]*rhouz[4]-0.31943828249997*Bz[4]*rhoux[4]-0.5*Bz[0]*rhoux[4]+0.31943828249997*Ey[4]*rho[4]+0.5*Ey[0]*rho[4]+0.5*rho[0]*Ey[4]-0.5*rhoux[0]*Bz[4]+0.5*rhouz[0]*Bx[4]+0.4472135954999579*Bx[3]*rhouz[3]-0.4472135954999579*Bz[3]*rhoux[3]+0.4472135954999579*Ey[3]*rho[3]+0.4472135954999579*Bx[1]*rhouz[1]-0.4472135954999579*Bz[1]*rhoux[1]+0.4472135954999579*Ey[1]*rho[1]; 
  outrhouy[5] += 0.31943828249997*Bx[7]*rhouz[7]+0.5000000000000001*Bx[1]*rhouz[7]-0.31943828249997*Bz[7]*rhoux[7]-0.5000000000000001*Bz[1]*rhoux[7]+0.31943828249997*Ey[7]*rho[7]+0.5000000000000001*Ey[1]*rho[7]+0.5000000000000001*rho[1]*Ey[7]-0.5000000000000001*rhoux[1]*Bz[7]+0.5000000000000001*rhouz[1]*Bx[7]+0.4472135954999579*Bx[6]*rhouz[6]-0.4472135954999579*Bz[6]*rhoux[6]+0.4472135954999579*Ey[6]*rho[6]+0.31943828249997*Bx[5]*rhouz[5]+0.5*Bx[0]*rhouz[5]-0.31943828249997*Bz[5]*rhoux[5]-0.5*Bz[0]*rhoux[5]+0.31943828249997*Ey[5]*rho[5]+0.5*Ey[0]*rho[5]+0.5*rho[0]*Ey[5]-0.5*rhoux[0]*Bz[5]+0.5*rhouz[0]*Bx[5]+0.4472135954999579*Bx[3]*rhouz[3]-0.4472135954999579*Bz[3]*rhoux[3]+0.4472135954999579*Ey[3]*rho[3]+0.4472135954999579*Bx[2]*rhouz[2]-0.4472135954999579*Bz[2]*rhoux[2]+0.4472135954999579*Ey[2]*rho[2]; 
  outrhouy[6] += 0.4*Bx[3]*rhouz[7]-0.4*Bz[3]*rhoux[7]+0.4*Ey[3]*rho[7]+0.4*rho[3]*Ey[7]-0.4*rhoux[3]*Bz[7]+0.4*rhouz[3]*Bx[7]+0.4472135954999579*Bx[5]*rhouz[6]+0.31943828249997*Bx[4]*rhouz[6]+0.5*Bx[0]*rhouz[6]-0.4472135954999579*Bz[5]*rhoux[6]-0.31943828249997*Bz[4]*rhoux[6]-0.5*Bz[0]*rhoux[6]+0.4472135954999579*Ey[5]*rho[6]+0.31943828249997*Ey[4]*rho[6]+0.5*Ey[0]*rho[6]+0.4472135954999579*rho[5]*Ey[6]+0.31943828249997*rho[4]*Ey[6]+0.5*rho[0]*Ey[6]-0.4472135954999579*rhoux[5]*Bz[6]-0.31943828249997*rhoux[4]*Bz[6]-0.5*rhoux[0]*Bz[6]+0.4472135954999579*rhouz[5]*Bx[6]+0.31943828249997*rhouz[4]*Bx[6]+0.5*rhouz[0]*Bx[6]+0.5000000000000001*Bx[2]*rhouz[4]-0.5000000000000001*Bz[2]*rhoux[4]+0.5000000000000001*Ey[2]*rho[4]+0.5000000000000001*rho[2]*Ey[4]-0.5000000000000001*rhoux[2]*Bz[4]+0.5000000000000001*rhouz[2]*Bx[4]+0.447213595499958*Bx[1]*rhouz[3]-0.447213595499958*Bz[1]*rhoux[3]+0.447213595499958*Ey[1]*rho[3]+0.447213595499958*rho[1]*Ey[3]-0.447213595499958*rhoux[1]*Bz[3]+0.447213595499958*rhouz[1]*Bx[3]; 
  outrhouy[7] += 0.31943828249997*Bx[5]*rhouz[7]+0.4472135954999579*Bx[4]*rhouz[7]+0.5*Bx[0]*rhouz[7]-0.31943828249997*Bz[5]*rhoux[7]-0.4472135954999579*Bz[4]*rhoux[7]-0.5*Bz[0]*rhoux[7]+0.31943828249997*Ey[5]*rho[7]+0.4472135954999579*Ey[4]*rho[7]+0.5*Ey[0]*rho[7]+0.31943828249997*rho[5]*Ey[7]+0.4472135954999579*rho[4]*Ey[7]+0.5*rho[0]*Ey[7]-0.31943828249997*rhoux[5]*Bz[7]-0.4472135954999579*rhoux[4]*Bz[7]-0.5*rhoux[0]*Bz[7]+0.31943828249997*rhouz[5]*Bx[7]+0.4472135954999579*rhouz[4]*Bx[7]+0.5*rhouz[0]*Bx[7]+0.4*Bx[3]*rhouz[6]-0.4*Bz[3]*rhoux[6]+0.4*Ey[3]*rho[6]+0.4*rho[3]*Ey[6]-0.4*rhoux[3]*Bz[6]+0.4*rhouz[3]*Bx[6]+0.5000000000000001*Bx[1]*rhouz[5]-0.5000000000000001*Bz[1]*rhoux[5]+0.5000000000000001*Ey[1]*rho[5]+0.5000000000000001*rho[1]*Ey[5]-0.5000000000000001*rhoux[1]*Bz[5]+0.5000000000000001*rhouz[1]*Bx[5]+0.447213595499958*Bx[2]*rhouz[3]-0.447213595499958*Bz[2]*rhoux[3]+0.447213595499958*Ey[2]*rho[3]+0.447213595499958*rho[2]*Ey[3]-0.447213595499958*rhoux[2]*Bz[3]+0.447213595499958*rhouz[2]*Bx[3]; 

  outrhouz[0] += (-0.5*Bx[7]*rhouy[7])+0.5*By[7]*rhoux[7]+0.5*Ez[7]*rho[7]-0.5*Bx[6]*rhouy[6]+0.5*By[6]*rhoux[6]+0.5*Ez[6]*rho[6]-0.5*Bx[5]*rhouy[5]+0.5*By[5]*rhoux[5]+0.5*Ez[5]*rho[5]-0.5*Bx[4]*rhouy[4]+0.5*By[4]*rhoux[4]+0.5*Ez[4]*rho[4]-0.5*Bx[3]*rhouy[3]+0.5*By[3]*rhoux[3]+0.5*Ez[3]*rho[3]-0.5*Bx[2]*rhouy[2]+0.5*By[2]*rhoux[2]+0.5*Ez[2]*rho[2]-0.5*Bx[1]*rhouy[1]+0.5*By[1]*rhoux[1]+0.5*Ez[1]*rho[1]-0.5*Bx[0]*rhouy[0]+0.5*By[0]*rhoux[0]+0.5*Ez[0]*rho[0]; 
  outrhouz[1] += (-0.5000000000000001*Bx[5]*rhouy[7])+0.5000000000000001*By[5]*rhoux[7]+0.5000000000000001*Ez[5]*rho[7]+0.5000000000000001*rho[5]*Ez[7]+0.5000000000000001*rhoux[5]*By[7]-0.5000000000000001*rhouy[5]*Bx[7]-0.447213595499958*Bx[3]*rhouy[6]+0.447213595499958*By[3]*rhoux[6]+0.447213595499958*Ez[3]*rho[6]+0.447213595499958*rho[3]*Ez[6]+0.447213595499958*rhoux[3]*By[6]-0.447213595499958*rhouy[3]*Bx[6]-0.4472135954999579*Bx[1]*rhouy[4]+0.4472135954999579*By[1]*rhoux[4]+0.4472135954999579*Ez[1]*rho[4]+0.4472135954999579*rho[1]*Ez[4]+0.4472135954999579*rhoux[1]*By[4]-0.4472135954999579*rhouy[1]*Bx[4]-0.5*Bx[2]*rhouy[3]+0.5*By[2]*rhoux[3]+0.5*Ez[2]*rho[3]+0.5*rho[2]*Ez[3]+0.5*rhoux[2]*By[3]-0.5*rhouy[2]*Bx[3]-0.5*Bx[0]*rhouy[1]+0.5*By[0]*rhoux[1]+0.5*Ez[0]*rho[1]+0.5*rho[0]*Ez[1]+0.5*rhoux[0]*By[1]-0.5*rhouy[0]*Bx[1]; 
  outrhouz[2] += (-0.447213595499958*Bx[3]*rhouy[7])+0.447213595499958*By[3]*rhoux[7]+0.447213595499958*Ez[3]*rho[7]+0.447213595499958*rho[3]*Ez[7]+0.447213595499958*rhoux[3]*By[7]-0.447213595499958*rhouy[3]*Bx[7]-0.5000000000000001*Bx[4]*rhouy[6]+0.5000000000000001*By[4]*rhoux[6]+0.5000000000000001*Ez[4]*rho[6]+0.5000000000000001*rho[4]*Ez[6]+0.5000000000000001*rhoux[4]*By[6]-0.5000000000000001*rhouy[4]*Bx[6]-0.4472135954999579*Bx[2]*rhouy[5]+0.4472135954999579*By[2]*rhoux[5]+0.4472135954999579*Ez[2]*rho[5]+0.4472135954999579*rho[2]*Ez[5]+0.4472135954999579*rhoux[2]*By[5]-0.4472135954999579*rhouy[2]*Bx[5]-0.5*Bx[1]*rhouy[3]+0.5*By[1]*rhoux[3]+0.5*Ez[1]*rho[3]+0.5*rho[1]*Ez[3]+0.5*rhoux[1]*By[3]-0.5*rhouy[1]*Bx[3]-0.5*Bx[0]*rhouy[2]+0.5*By[0]*rhoux[2]+0.5*Ez[0]*rho[2]+0.5*rho[0]*Ez[2]+0.5*rhoux[0]*By[2]-0.5*rhouy[0]*Bx[2]; 
  outrhouz[3] += (-0.4*Bx[6]*rhouy[7])-0.447213595499958*Bx[2]*rhouy[7]+0.4*By[6]*rhoux[7]+0.447213595499958*By[2]*rhoux[7]+0.4*Ez[6]*rho[7]+0.447213595499958*Ez[2]*rho[7]+0.4*rho[6]*Ez[7]+0.447213595499958*rho[2]*Ez[7]+0.4*rhoux[6]*By[7]+0.447213595499958*rhoux[2]*By[7]-0.4*rhouy[6]*Bx[7]-0.447213595499958*rhouy[2]*Bx[7]-0.447213595499958*Bx[1]*rhouy[6]+0.447213595499958*By[1]*rhoux[6]+0.447213595499958*Ez[1]*rho[6]+0.447213595499958*rho[1]*Ez[6]+0.447213595499958*rhoux[1]*By[6]-0.447213595499958*rhouy[1]*Bx[6]-0.4472135954999579*Bx[3]*rhouy[5]+0.4472135954999579*By[3]*rhoux[5]+0.4472135954999579*Ez[3]*rho[5]+0.4472135954999579*rho[3]*Ez[5]+0.4472135954999579*rhoux[3]*By[5]-0.4472135954999579*rhouy[3]*Bx[5]-0.4472135954999579*Bx[3]*rhouy[4]+0.4472135954999579*By[3]*rhoux[4]+0.4472135954999579*Ez[3]*rho[4]+0.4472135954999579*rho[3]*Ez[4]+0.4472135954999579*rhoux[3]*By[4]-0.4472135954999579*rhouy[3]*Bx[4]-0.5*Bx[0]*rhouy[3]+0.5*By[0]*rhoux[3]+0.5*Ez[0]*rho[3]+0.5*rho[0]*Ez[3]+0.5*rhoux[0]*By[3]-0.5*rhouy[0]*Bx[3]-0.5*Bx[1]*rhouy[2]+0.5*By[1]*rhoux[2]+0.5*Ez[1]*rho[2]+0.5*rho[1]*Ez[2]+0.5*rhoux[1]*By[2]-0.5*rhouy[1]*Bx[2]; 
  outrhouz[4] += (-0.4472135954999579*Bx[7]*rhouy[7])+0.4472135954999579*By[7]*rhoux[7]+0.4472135954999579*Ez[7]*rho[7]-0.31943828249997*Bx[6]*rhouy[6]-0.5000000000000001*Bx[2]*rhouy[6]+0.31943828249997*By[6]*rhoux[6]+0.5000000000000001*By[2]*rhoux[6]+0.31943828249997*Ez[6]*rho[6]+0.5000000000000001*Ez[2]*rho[6]+0.5000000000000001*rho[2]*Ez[6]+0.5000000000000001*rhoux[2]*By[6]-0.5000000000000001*rhouy[2]*Bx[6]-0.31943828249997*Bx[4]*rhouy[4]-0.5*Bx[0]*rhouy[4]+0.31943828249997*By[4]*rhoux[4]+0.5*By[0]*rhoux[4]+0.31943828249997*Ez[4]*rho[4]+0.5*Ez[0]*rho[4]+0.5*rho[0]*Ez[4]+0.5*rhoux[0]*By[4]-0.5*rhouy[0]*Bx[4]-0.4472135954999579*Bx[3]*rhouy[3]+0.4472135954999579*By[3]*rhoux[3]+0.4472135954999579*Ez[3]*rho[3]-0.4472135954999579*Bx[1]*rhouy[1]+0.4472135954999579*By[1]*rhoux[1]+0.4472135954999579*Ez[1]*rho[1]; 
  outrhouz[5] += (-0.31943828249997*Bx[7]*rhouy[7])-0.5000000000000001*Bx[1]*rhouy[7]+0.31943828249997*By[7]*rhoux[7]+0.5000000000000001*By[1]*rhoux[7]+0.31943828249997*Ez[7]*rho[7]+0.5000000000000001*Ez[1]*rho[7]+0.5000000000000001*rho[1]*Ez[7]+0.5000000000000001*rhoux[1]*By[7]-0.5000000000000001*rhouy[1]*Bx[7]-0.4472135954999579*Bx[6]*rhouy[6]+0.4472135954999579*By[6]*rhoux[6]+0.4472135954999579*Ez[6]*rho[6]-0.31943828249997*Bx[5]*rhouy[5]-0.5*Bx[0]*rhouy[5]+0.31943828249997*By[5]*rhoux[5]+0.5*By[0]*rhoux[5]+0.31943828249997*Ez[5]*rho[5]+0.5*Ez[0]*rho[5]+0.5*rho[0]*Ez[5]+0.5*rhoux[0]*By[5]-0.5*rhouy[0]*Bx[5]-0.4472135954999579*Bx[3]*rhouy[3]+0.4472135954999579*By[3]*rhoux[3]+0.4472135954999579*Ez[3]*rho[3]-0.4472135954999579*Bx[2]*rhouy[2]+0.4472135954999579*By[2]*rhoux[2]+0.4472135954999579*Ez[2]*rho[2]; 
  outrhouz[6] += (-0.4*Bx[3]*rhouy[7])+0.4*By[3]*rhoux[7]+0.4*Ez[3]*rho[7]+0.4*rho[3]*Ez[7]+0.4*rhoux[3]*By[7]-0.4*rhouy[3]*Bx[7]-0.4472135954999579*Bx[5]*rhouy[6]-0.31943828249997*Bx[4]*rhouy[6]-0.5*Bx[0]*rhouy[6]+0.4472135954999579*By[5]*rhoux[6]+0.31943828249997*By[4]*rhoux[6]+0.5*By[0]*rhoux[6]+0.4472135954999579*Ez[5]*rho[6]+0.31943828249997*Ez[4]*rho[6]+0.5*Ez[0]*rho[6]+0.4472135954999579*rho[5]*Ez[6]+0.31943828249997*rho[4]*Ez[6]+0.5*rho[0]*Ez[6]+0.4472135954999579*rhoux[5]*By[6]+0.31943828249997*rhoux[4]*By[6]+0.5*rhoux[0]*By[6]-0.4472135954999579*rhouy[5]*Bx[6]-0.31943828249997*rhouy[4]*Bx[6]-0.5*rhouy[0]*Bx[6]-0.5000000000000001*Bx[2]*rhouy[4]+0.5000000000000001*By[2]*rhoux[4]+0.5000000000000001*Ez[2]*rho[4]+0.5000000000000001*rho[2]*Ez[4]+0.5000000000000001*rhoux[2]*By[4]-0.5000000000000001*rhouy[2]*Bx[4]-0.447213595499958*Bx[1]*rhouy[3]+0.447213595499958*By[1]*rhoux[3]+0.447213595499958*Ez[1]*rho[3]+0.447213595499958*rho[1]*Ez[3]+0.447213595499958*rhoux[1]*By[3]-0.447213595499958*rhouy[1]*Bx[3]; 
  outrhouz[7] += (-0.31943828249997*Bx[5]*rhouy[7])-0.4472135954999579*Bx[4]*rhouy[7]-0.5*Bx[0]*rhouy[7]+0.31943828249997*By[5]*rhoux[7]+0.4472135954999579*By[4]*rhoux[7]+0.5*By[0]*rhoux[7]+0.31943828249997*Ez[5]*rho[7]+0.4472135954999579*Ez[4]*rho[7]+0.5*Ez[0]*rho[7]+0.31943828249997*rho[5]*Ez[7]+0.4472135954999579*rho[4]*Ez[7]+0.5*rho[0]*Ez[7]+0.31943828249997*rhoux[5]*By[7]+0.4472135954999579*rhoux[4]*By[7]+0.5*rhoux[0]*By[7]-0.31943828249997*rhouy[5]*Bx[7]-0.4472135954999579*rhouy[4]*Bx[7]-0.5*rhouy[0]*Bx[7]-0.4*Bx[3]*rhouy[6]+0.4*By[3]*rhoux[6]+0.4*Ez[3]*rho[6]+0.4*rho[3]*Ez[6]+0.4*rhoux[3]*By[6]-0.4*rhouy[3]*Bx[6]-0.5000000000000001*Bx[1]*rhouy[5]+0.5000000000000001*By[1]*rhoux[5]+0.5000000000000001*Ez[1]*rho[5]+0.5000000000000001*rho[1]*Ez[5]+0.5000000000000001*rhoux[1]*By[5]-0.5000000000000001*rhouy[1]*Bx[5]-0.447213595499958*Bx[2]*rhouy[3]+0.447213595499958*By[2]*rhoux[3]+0.447213595499958*Ez[2]*rho[3]+0.447213595499958*rho[2]*Ez[3]+0.447213595499958*rhoux[2]*By[3]-0.447213595499958*rhouy[2]*Bx[3]; 

  outE_perp[0] += 0.5*Ez[7]*rhou_perp_z[7]+0.5*Ey[7]*rhou_perp_y[7]+0.5*Ex[7]*rhou_perp_x[7]+0.5*nu_vth_sq[7]*rho[7]-0.5*nu[7]*p_perp[7]+0.5*Ez[6]*rhou_perp_z[6]+0.5*Ey[6]*rhou_perp_y[6]+0.5*Ex[6]*rhou_perp_x[6]+0.5*nu_vth_sq[6]*rho[6]-0.5*nu[6]*p_perp[6]+0.5*Ez[5]*rhou_perp_z[5]+0.5*Ey[5]*rhou_perp_y[5]+0.5*Ex[5]*rhou_perp_x[5]+0.5*nu_vth_sq[5]*rho[5]-0.5*nu[5]*p_perp[5]+0.5*Ez[4]*rhou_perp_z[4]+0.5*Ey[4]*rhou_perp_y[4]+0.5*Ex[4]*rhou_perp_x[4]+0.5*nu_vth_sq[4]*rho[4]-0.5*nu[4]*p_perp[4]+0.5*Ez[3]*rhou_perp_z[3]+0.5*Ey[3]*rhou_perp_y[3]+0.5*Ex[3]*rhou_perp_x[3]+0.5*nu_vth_sq[3]*rho[3]-0.5*nu[3]*p_perp[3]+0.5*Ez[2]*rhou_perp_z[2]+0.5*Ey[2]*rhou_perp_y[2]+0.5*Ex[2]*rhou_perp_x[2]+0.5*nu_vth_sq[2]*rho[2]-0.5*nu[2]*p_perp[2]+0.5*Ez[1]*rhou_perp_z[1]+0.5*Ey[1]*rhou_perp_y[1]+0.5*Ex[1]*rhou_perp_x[1]+0.5*nu_vth_sq[1]*rho[1]-0.5*nu[1]*p_perp[1]+0.5*Ez[0]*rhou_perp_z[0]+0.5*Ey[0]*rhou_perp_y[0]+0.5*Ex[0]*rhou_perp_x[0]+0.5*nu_vth_sq[0]*rho[0]-0.5*nu[0]*p_perp[0]; 
  outE_perp[1] += 0.5000000000000001*Ez[5]*rhou_perp_z[7]+0.5000000000000001*Ey[5]*rhou_perp_y[7]+0.5000000000000001*Ex[5]*rhou_perp_x[7]+0.5000000000000001*nu_vth_sq[5]*rho[7]-0.5000000000000001*nu[5]*p_perp[7]+0.5000000000000001*rho[5]*nu_vth_sq[7]-0.5000000000000001*p_perp[5]*nu[7]+0.5000000000000001*rhou_perp_z[5]*Ez[7]+0.5000000000000001*rhou_perp_y[5]*Ey[7]+0.5000000000000001*rhou_perp_x[5]*Ex[7]+0.447213595499958*Ez[3]*rhou_perp_z[6]+0.447213595499958*Ey[3]*rhou_perp_y[6]+0.447213595499958*Ex[3]*rhou_perp_x[6]+0.447213595499958*nu_vth_sq[3]*rho[6]-0.447213595499958*nu[3]*p_perp[6]+0.447213595499958*rho[3]*nu_vth_sq[6]-0.447213595499958*p_perp[3]*nu[6]+0.447213595499958*rhou_perp_z[3]*Ez[6]+0.447213595499958*rhou_perp_y[3]*Ey[6]+0.447213595499958*rhou_perp_x[3]*Ex[6]+0.4472135954999579*Ez[1]*rhou_perp_z[4]+0.4472135954999579*Ey[1]*rhou_perp_y[4]+0.4472135954999579*Ex[1]*rhou_perp_x[4]+0.4472135954999579*nu_vth_sq[1]*rho[4]-0.4472135954999579*nu[1]*p_perp[4]+0.4472135954999579*rho[1]*nu_vth_sq[4]-0.4472135954999579*p_perp[1]*nu[4]+0.4472135954999579*rhou_perp_z[1]*Ez[4]+0.4472135954999579*rhou_perp_y[1]*Ey[4]+0.4472135954999579*rhou_perp_x[1]*Ex[4]+0.5*Ez[2]*rhou_perp_z[3]+0.5*Ey[2]*rhou_perp_y[3]+0.5*Ex[2]*rhou_perp_x[3]+0.5*nu_vth_sq[2]*rho[3]-0.5*nu[2]*p_perp[3]+0.5*rho[2]*nu_vth_sq[3]-0.5*p_perp[2]*nu[3]+0.5*rhou_perp_z[2]*Ez[3]+0.5*rhou_perp_y[2]*Ey[3]+0.5*rhou_perp_x[2]*Ex[3]+0.5*Ez[0]*rhou_perp_z[1]+0.5*Ey[0]*rhou_perp_y[1]+0.5*Ex[0]*rhou_perp_x[1]+0.5*nu_vth_sq[0]*rho[1]-0.5*nu[0]*p_perp[1]+0.5*rho[0]*nu_vth_sq[1]-0.5*p_perp[0]*nu[1]+0.5*rhou_perp_z[0]*Ez[1]+0.5*rhou_perp_y[0]*Ey[1]+0.5*rhou_perp_x[0]*Ex[1]; 
  outE_perp[2] += 0.447213595499958*Ez[3]*rhou_perp_z[7]+0.447213595499958*Ey[3]*rhou_perp_y[7]+0.447213595499958*Ex[3]*rhou_perp_x[7]+0.447213595499958*nu_vth_sq[3]*rho[7]-0.447213595499958*nu[3]*p_perp[7]+0.447213595499958*rho[3]*nu_vth_sq[7]-0.447213595499958*p_perp[3]*nu[7]+0.447213595499958*rhou_perp_z[3]*Ez[7]+0.447213595499958*rhou_perp_y[3]*Ey[7]+0.447213595499958*rhou_perp_x[3]*Ex[7]+0.5000000000000001*Ez[4]*rhou_perp_z[6]+0.5000000000000001*Ey[4]*rhou_perp_y[6]+0.5000000000000001*Ex[4]*rhou_perp_x[6]+0.5000000000000001*nu_vth_sq[4]*rho[6]-0.5000000000000001*nu[4]*p_perp[6]+0.5000000000000001*rho[4]*nu_vth_sq[6]-0.5000000000000001*p_perp[4]*nu[6]+0.5000000000000001*rhou_perp_z[4]*Ez[6]+0.5000000000000001*rhou_perp_y[4]*Ey[6]+0.5000000000000001*rhou_perp_x[4]*Ex[6]+0.4472135954999579*Ez[2]*rhou_perp_z[5]+0.4472135954999579*Ey[2]*rhou_perp_y[5]+0.4472135954999579*Ex[2]*rhou_perp_x[5]+0.4472135954999579*nu_vth_sq[2]*rho[5]-0.4472135954999579*nu[2]*p_perp[5]+0.4472135954999579*rho[2]*nu_vth_sq[5]-0.4472135954999579*p_perp[2]*nu[5]+0.4472135954999579*rhou_perp_z[2]*Ez[5]+0.4472135954999579*rhou_perp_y[2]*Ey[5]+0.4472135954999579*rhou_perp_x[2]*Ex[5]+0.5*Ez[1]*rhou_perp_z[3]+0.5*Ey[1]*rhou_perp_y[3]+0.5*Ex[1]*rhou_perp_x[3]+0.5*nu_vth_sq[1]*rho[3]-0.5*nu[1]*p_perp[3]+0.5*rho[1]*nu_vth_sq[3]-0.5*p_perp[1]*nu[3]+0.5*rhou_perp_z[1]*Ez[3]+0.5*rhou_perp_y[1]*Ey[3]+0.5*rhou_perp_x[1]*Ex[3]+0.5*Ez[0]*rhou_perp_z[2]+0.5*Ey[0]*rhou_perp_y[2]+0.5*Ex[0]*rhou_perp_x[2]+0.5*nu_vth_sq[0]*rho[2]-0.5*nu[0]*p_perp[2]+0.5*rho[0]*nu_vth_sq[2]-0.5*p_perp[0]*nu[2]+0.5*rhou_perp_z[0]*Ez[2]+0.5*rhou_perp_y[0]*Ey[2]+0.5*rhou_perp_x[0]*Ex[2]; 
  outE_perp[3] += 0.4*Ez[6]*rhou_perp_z[7]+0.447213595499958*Ez[2]*rhou_perp_z[7]+0.4*Ey[6]*rhou_perp_y[7]+0.447213595499958*Ey[2]*rhou_perp_y[7]+0.4*Ex[6]*rhou_perp_x[7]+0.447213595499958*Ex[2]*rhou_perp_x[7]+0.4*nu_vth_sq[6]*rho[7]+0.447213595499958*nu_vth_sq[2]*rho[7]-0.4*nu[6]*p_perp[7]-0.447213595499958*nu[2]*p_perp[7]+0.4*rho[6]*nu_vth_sq[7]+0.447213595499958*rho[2]*nu_vth_sq[7]-0.4*p_perp[6]*nu[7]-0.447213595499958*p_perp[2]*nu[7]+0.4*rhou_perp_z[6]*Ez[7]+0.447213595499958*rhou_perp_z[2]*Ez[7]+0.4*rhou_perp_y[6]*Ey[7]+0.447213595499958*rhou_perp_y[2]*Ey[7]+0.4*rhou_perp_x[6]*Ex[7]+0.447213595499958*rhou_perp_x[2]*Ex[7]+0.447213595499958*Ez[1]*rhou_perp_z[6]+0.447213595499958*Ey[1]*rhou_perp_y[6]+0.447213595499958*Ex[1]*rhou_perp_x[6]+0.447213595499958*nu_vth_sq[1]*rho[6]-0.447213595499958*nu[1]*p_perp[6]+0.447213595499958*rho[1]*nu_vth_sq[6]-0.447213595499958*p_perp[1]*nu[6]+0.447213595499958*rhou_perp_z[1]*Ez[6]+0.447213595499958*rhou_perp_y[1]*Ey[6]+0.447213595499958*rhou_perp_x[1]*Ex[6]+0.4472135954999579*Ez[3]*rhou_perp_z[5]+0.4472135954999579*Ey[3]*rhou_perp_y[5]+0.4472135954999579*Ex[3]*rhou_perp_x[5]+0.4472135954999579*nu_vth_sq[3]*rho[5]-0.4472135954999579*nu[3]*p_perp[5]+0.4472135954999579*rho[3]*nu_vth_sq[5]-0.4472135954999579*p_perp[3]*nu[5]+0.4472135954999579*rhou_perp_z[3]*Ez[5]+0.4472135954999579*rhou_perp_y[3]*Ey[5]+0.4472135954999579*rhou_perp_x[3]*Ex[5]+0.4472135954999579*Ez[3]*rhou_perp_z[4]+0.4472135954999579*Ey[3]*rhou_perp_y[4]+0.4472135954999579*Ex[3]*rhou_perp_x[4]+0.4472135954999579*nu_vth_sq[3]*rho[4]-0.4472135954999579*nu[3]*p_perp[4]+0.4472135954999579*rho[3]*nu_vth_sq[4]-0.4472135954999579*p_perp[3]*nu[4]+0.4472135954999579*rhou_perp_z[3]*Ez[4]+0.4472135954999579*rhou_perp_y[3]*Ey[4]+0.4472135954999579*rhou_perp_x[3]*Ex[4]+0.5*Ez[0]*rhou_perp_z[3]+0.5*Ey[0]*rhou_perp_y[3]+0.5*Ex[0]*rhou_perp_x[3]+0.5*nu_vth_sq[0]*rho[3]-0.5*nu[0]*p_perp[3]+0.5*rho[0]*nu_vth_sq[3]-0.5*p_perp[0]*nu[3]+0.5*rhou_perp_z[0]*Ez[3]+0.5*rhou_perp_y[0]*Ey[3]+0.5*rhou_perp_x[0]*Ex[3]+0.5*Ez[1]*rhou_perp_z[2]+0.5*Ey[1]*rhou_perp_y[2]+0.5*Ex[1]*rhou_perp_x[2]+0.5*nu_vth_sq[1]*rho[2]-0.5*nu[1]*p_perp[2]+0.5*rho[1]*nu_vth_sq[2]-0.5*p_perp[1]*nu[2]+0.5*rhou_perp_z[1]*Ez[2]+0.5*rhou_perp_y[1]*Ey[2]+0.5*rhou_perp_x[1]*Ex[2]; 
  outE_perp[4] += 0.4472135954999579*Ez[7]*rhou_perp_z[7]+0.4472135954999579*Ey[7]*rhou_perp_y[7]+0.4472135954999579*Ex[7]*rhou_perp_x[7]+0.4472135954999579*nu_vth_sq[7]*rho[7]-0.4472135954999579*nu[7]*p_perp[7]+0.31943828249997*Ez[6]*rhou_perp_z[6]+0.5000000000000001*Ez[2]*rhou_perp_z[6]+0.31943828249997*Ey[6]*rhou_perp_y[6]+0.5000000000000001*Ey[2]*rhou_perp_y[6]+0.31943828249997*Ex[6]*rhou_perp_x[6]+0.5000000000000001*Ex[2]*rhou_perp_x[6]+0.31943828249997*nu_vth_sq[6]*rho[6]+0.5000000000000001*nu_vth_sq[2]*rho[6]-0.31943828249997*nu[6]*p_perp[6]-0.5000000000000001*nu[2]*p_perp[6]+0.5000000000000001*rho[2]*nu_vth_sq[6]-0.5000000000000001*p_perp[2]*nu[6]+0.5000000000000001*rhou_perp_z[2]*Ez[6]+0.5000000000000001*rhou_perp_y[2]*Ey[6]+0.5000000000000001*rhou_perp_x[2]*Ex[6]+0.31943828249997*Ez[4]*rhou_perp_z[4]+0.5*Ez[0]*rhou_perp_z[4]+0.31943828249997*Ey[4]*rhou_perp_y[4]+0.5*Ey[0]*rhou_perp_y[4]+0.31943828249997*Ex[4]*rhou_perp_x[4]+0.5*Ex[0]*rhou_perp_x[4]+0.31943828249997*nu_vth_sq[4]*rho[4]+0.5*nu_vth_sq[0]*rho[4]-0.31943828249997*nu[4]*p_perp[4]-0.5*nu[0]*p_perp[4]+0.5*rho[0]*nu_vth_sq[4]-0.5*p_perp[0]*nu[4]+0.5*rhou_perp_z[0]*Ez[4]+0.5*rhou_perp_y[0]*Ey[4]+0.5*rhou_perp_x[0]*Ex[4]+0.4472135954999579*Ez[3]*rhou_perp_z[3]+0.4472135954999579*Ey[3]*rhou_perp_y[3]+0.4472135954999579*Ex[3]*rhou_perp_x[3]+0.4472135954999579*nu_vth_sq[3]*rho[3]-0.4472135954999579*nu[3]*p_perp[3]+0.4472135954999579*Ez[1]*rhou_perp_z[1]+0.4472135954999579*Ey[1]*rhou_perp_y[1]+0.4472135954999579*Ex[1]*rhou_perp_x[1]+0.4472135954999579*nu_vth_sq[1]*rho[1]-0.4472135954999579*nu[1]*p_perp[1]; 
  outE_perp[5] += 0.31943828249997*Ez[7]*rhou_perp_z[7]+0.5000000000000001*Ez[1]*rhou_perp_z[7]+0.31943828249997*Ey[7]*rhou_perp_y[7]+0.5000000000000001*Ey[1]*rhou_perp_y[7]+0.31943828249997*Ex[7]*rhou_perp_x[7]+0.5000000000000001*Ex[1]*rhou_perp_x[7]+0.31943828249997*nu_vth_sq[7]*rho[7]+0.5000000000000001*nu_vth_sq[1]*rho[7]-0.31943828249997*nu[7]*p_perp[7]-0.5000000000000001*nu[1]*p_perp[7]+0.5000000000000001*rho[1]*nu_vth_sq[7]-0.5000000000000001*p_perp[1]*nu[7]+0.5000000000000001*rhou_perp_z[1]*Ez[7]+0.5000000000000001*rhou_perp_y[1]*Ey[7]+0.5000000000000001*rhou_perp_x[1]*Ex[7]+0.4472135954999579*Ez[6]*rhou_perp_z[6]+0.4472135954999579*Ey[6]*rhou_perp_y[6]+0.4472135954999579*Ex[6]*rhou_perp_x[6]+0.4472135954999579*nu_vth_sq[6]*rho[6]-0.4472135954999579*nu[6]*p_perp[6]+0.31943828249997*Ez[5]*rhou_perp_z[5]+0.5*Ez[0]*rhou_perp_z[5]+0.31943828249997*Ey[5]*rhou_perp_y[5]+0.5*Ey[0]*rhou_perp_y[5]+0.31943828249997*Ex[5]*rhou_perp_x[5]+0.5*Ex[0]*rhou_perp_x[5]+0.31943828249997*nu_vth_sq[5]*rho[5]+0.5*nu_vth_sq[0]*rho[5]-0.31943828249997*nu[5]*p_perp[5]-0.5*nu[0]*p_perp[5]+0.5*rho[0]*nu_vth_sq[5]-0.5*p_perp[0]*nu[5]+0.5*rhou_perp_z[0]*Ez[5]+0.5*rhou_perp_y[0]*Ey[5]+0.5*rhou_perp_x[0]*Ex[5]+0.4472135954999579*Ez[3]*rhou_perp_z[3]+0.4472135954999579*Ey[3]*rhou_perp_y[3]+0.4472135954999579*Ex[3]*rhou_perp_x[3]+0.4472135954999579*nu_vth_sq[3]*rho[3]-0.4472135954999579*nu[3]*p_perp[3]+0.4472135954999579*Ez[2]*rhou_perp_z[2]+0.4472135954999579*Ey[2]*rhou_perp_y[2]+0.4472135954999579*Ex[2]*rhou_perp_x[2]+0.4472135954999579*nu_vth_sq[2]*rho[2]-0.4472135954999579*nu[2]*p_perp[2]; 
  outE_perp[6] += 0.4*Ez[3]*rhou_perp_z[7]+0.4*Ey[3]*rhou_perp_y[7]+0.4*Ex[3]*rhou_perp_x[7]+0.4*nu_vth_sq[3]*rho[7]-0.4*nu[3]*p_perp[7]+0.4*rho[3]*nu_vth_sq[7]-0.4*p_perp[3]*nu[7]+0.4*rhou_perp_z[3]*Ez[7]+0.4*rhou_perp_y[3]*Ey[7]+0.4*rhou_perp_x[3]*Ex[7]+0.4472135954999579*Ez[5]*rhou_perp_z[6]+0.31943828249997*Ez[4]*rhou_perp_z[6]+0.5*Ez[0]*rhou_perp_z[6]+0.4472135954999579*Ey[5]*rhou_perp_y[6]+0.31943828249997*Ey[4]*rhou_perp_y[6]+0.5*Ey[0]*rhou_perp_y[6]+0.4472135954999579*Ex[5]*rhou_perp_x[6]+0.31943828249997*Ex[4]*rhou_perp_x[6]+0.5*Ex[0]*rhou_perp_x[6]+0.4472135954999579*nu_vth_sq[5]*rho[6]+0.31943828249997*nu_vth_sq[4]*rho[6]+0.5*nu_vth_sq[0]*rho[6]-0.4472135954999579*nu[5]*p_perp[6]-0.31943828249997*nu[4]*p_perp[6]-0.5*nu[0]*p_perp[6]+0.4472135954999579*rho[5]*nu_vth_sq[6]+0.31943828249997*rho[4]*nu_vth_sq[6]+0.5*rho[0]*nu_vth_sq[6]-0.4472135954999579*p_perp[5]*nu[6]-0.31943828249997*p_perp[4]*nu[6]-0.5*p_perp[0]*nu[6]+0.4472135954999579*rhou_perp_z[5]*Ez[6]+0.31943828249997*rhou_perp_z[4]*Ez[6]+0.5*rhou_perp_z[0]*Ez[6]+0.4472135954999579*rhou_perp_y[5]*Ey[6]+0.31943828249997*rhou_perp_y[4]*Ey[6]+0.5*rhou_perp_y[0]*Ey[6]+0.4472135954999579*rhou_perp_x[5]*Ex[6]+0.31943828249997*rhou_perp_x[4]*Ex[6]+0.5*rhou_perp_x[0]*Ex[6]+0.5000000000000001*Ez[2]*rhou_perp_z[4]+0.5000000000000001*Ey[2]*rhou_perp_y[4]+0.5000000000000001*Ex[2]*rhou_perp_x[4]+0.5000000000000001*nu_vth_sq[2]*rho[4]-0.5000000000000001*nu[2]*p_perp[4]+0.5000000000000001*rho[2]*nu_vth_sq[4]-0.5000000000000001*p_perp[2]*nu[4]+0.5000000000000001*rhou_perp_z[2]*Ez[4]+0.5000000000000001*rhou_perp_y[2]*Ey[4]+0.5000000000000001*rhou_perp_x[2]*Ex[4]+0.447213595499958*Ez[1]*rhou_perp_z[3]+0.447213595499958*Ey[1]*rhou_perp_y[3]+0.447213595499958*Ex[1]*rhou_perp_x[3]+0.447213595499958*nu_vth_sq[1]*rho[3]-0.447213595499958*nu[1]*p_perp[3]+0.447213595499958*rho[1]*nu_vth_sq[3]-0.447213595499958*p_perp[1]*nu[3]+0.447213595499958*rhou_perp_z[1]*Ez[3]+0.447213595499958*rhou_perp_y[1]*Ey[3]+0.447213595499958*rhou_perp_x[1]*Ex[3]; 
  outE_perp[7] += 0.31943828249997*Ez[5]*rhou_perp_z[7]+0.4472135954999579*Ez[4]*rhou_perp_z[7]+0.5*Ez[0]*rhou_perp_z[7]+0.31943828249997*Ey[5]*rhou_perp_y[7]+0.4472135954999579*Ey[4]*rhou_perp_y[7]+0.5*Ey[0]*rhou_perp_y[7]+0.31943828249997*Ex[5]*rhou_perp_x[7]+0.4472135954999579*Ex[4]*rhou_perp_x[7]+0.5*Ex[0]*rhou_perp_x[7]+0.31943828249997*nu_vth_sq[5]*rho[7]+0.4472135954999579*nu_vth_sq[4]*rho[7]+0.5*nu_vth_sq[0]*rho[7]-0.31943828249997*nu[5]*p_perp[7]-0.4472135954999579*nu[4]*p_perp[7]-0.5*nu[0]*p_perp[7]+0.31943828249997*rho[5]*nu_vth_sq[7]+0.4472135954999579*rho[4]*nu_vth_sq[7]+0.5*rho[0]*nu_vth_sq[7]-0.31943828249997*p_perp[5]*nu[7]-0.4472135954999579*p_perp[4]*nu[7]-0.5*p_perp[0]*nu[7]+0.31943828249997*rhou_perp_z[5]*Ez[7]+0.4472135954999579*rhou_perp_z[4]*Ez[7]+0.5*rhou_perp_z[0]*Ez[7]+0.31943828249997*rhou_perp_y[5]*Ey[7]+0.4472135954999579*rhou_perp_y[4]*Ey[7]+0.5*rhou_perp_y[0]*Ey[7]+0.31943828249997*rhou_perp_x[5]*Ex[7]+0.4472135954999579*rhou_perp_x[4]*Ex[7]+0.5*rhou_perp_x[0]*Ex[7]+0.4*Ez[3]*rhou_perp_z[6]+0.4*Ey[3]*rhou_perp_y[6]+0.4*Ex[3]*rhou_perp_x[6]+0.4*nu_vth_sq[3]*rho[6]-0.4*nu[3]*p_perp[6]+0.4*rho[3]*nu_vth_sq[6]-0.4*p_perp[3]*nu[6]+0.4*rhou_perp_z[3]*Ez[6]+0.4*rhou_perp_y[3]*Ey[6]+0.4*rhou_perp_x[3]*Ex[6]+0.5000000000000001*Ez[1]*rhou_perp_z[5]+0.5000000000000001*Ey[1]*rhou_perp_y[5]+0.5000000000000001*Ex[1]*rhou_perp_x[5]+0.5000000000000001*nu_vth_sq[1]*rho[5]-0.5000000000000001*nu[1]*p_perp[5]+0.5000000000000001*rho[1]*nu_vth_sq[5]-0.5000000000000001*p_perp[1]*nu[5]+0.5000000000000001*rhou_perp_z[1]*Ez[5]+0.5000000000000001*rhou_perp_y[1]*Ey[5]+0.5000000000000001*rhou_perp_x[1]*Ex[5]+0.447213595499958*Ez[2]*rhou_perp_z[3]+0.447213595499958*Ey[2]*rhou_perp_y[3]+0.447213595499958*Ex[2]*rhou_perp_x[3]+0.447213595499958*nu_vth_sq[2]*rho[3]-0.447213595499958*nu[2]*p_perp[3]+0.447213595499958*rho[2]*nu_vth_sq[3]-0.447213595499958*p_perp[2]*nu[3]+0.447213595499958*rhou_perp_z[2]*Ez[3]+0.447213595499958*rhou_perp_y[2]*Ey[3]+0.447213595499958*rhou_perp_x[2]*Ex[3]; 

} 
