#include <gkyl_mat.h> 
#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void fluid_vars_u_copy_2x_ser_p1(int count, struct gkyl_nmat *x, 
    double* GKYL_RESTRICT u, double* GKYL_RESTRICT u_surf) 
{ 
  // count:  integer to indicate which matrix being fetched. 
  // x:      Input solution vector. 
  // u:      Output volume expansion of flow velocity: 
  //         [ux, uy, uz]. 
  // u_surf: Output surface expansion of flow velocity 
  //         [ux_xl, ux_xr, uy_xl, uy_xr, uz_xl, uz_xr, 
  //          ux_yl, ux_yr, uy_yl, uy_yr, uz_yl, uz_yr,  
  //          ux_zl, ux_zr, uy_zl, uy_zr, uz_zl, uz_zr]  
 
  struct gkyl_mat x_ux = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_uy = gkyl_nmat_get(x, count+1); 
  struct gkyl_mat x_uz = gkyl_nmat_get(x, count+2); 
  double *ux = &u[0]; 
  double *uy = &u[4]; 
  double *uz = &u[8]; 

  ux[0] = gkyl_mat_get(&x_ux,0,0); 
  uy[0] = gkyl_mat_get(&x_uy,0,0); 
  uz[0] = gkyl_mat_get(&x_uz,0,0); 
  ux[1] = gkyl_mat_get(&x_ux,1,0); 
  uy[1] = gkyl_mat_get(&x_uy,1,0); 
  uz[1] = gkyl_mat_get(&x_uz,1,0); 
  ux[2] = gkyl_mat_get(&x_ux,2,0); 
  uy[2] = gkyl_mat_get(&x_uy,2,0); 
  uz[2] = gkyl_mat_get(&x_uz,2,0); 
  ux[3] = gkyl_mat_get(&x_ux,3,0); 
  uy[3] = gkyl_mat_get(&x_uy,3,0); 
  uz[3] = gkyl_mat_get(&x_uz,3,0); 

  double *ux_xl = &u_surf[0]; 
  double *ux_xr = &u_surf[2]; 
  double *uy_xl = &u_surf[4]; 
  double *uy_xr = &u_surf[6]; 
  double *uz_xl = &u_surf[8]; 
  double *uz_xr = &u_surf[10]; 
 
  ux_xl[0] = 0.7071067811865475*ux[0]-1.224744871391589*ux[1]; 
  ux_xl[1] = 0.7071067811865475*ux[2]-1.224744871391589*ux[3]; 
  uy_xl[0] = 0.7071067811865475*uy[0]-1.224744871391589*uy[1]; 
  uy_xl[1] = 0.7071067811865475*uy[2]-1.224744871391589*uy[3]; 
  uz_xl[0] = 0.7071067811865475*uz[0]-1.224744871391589*uz[1]; 
  uz_xl[1] = 0.7071067811865475*uz[2]-1.224744871391589*uz[3]; 
 
  ux_xr[0] = 1.224744871391589*ux[1]+0.7071067811865475*ux[0]; 
  ux_xr[1] = 1.224744871391589*ux[3]+0.7071067811865475*ux[2]; 
  uy_xr[0] = 1.224744871391589*uy[1]+0.7071067811865475*uy[0]; 
  uy_xr[1] = 1.224744871391589*uy[3]+0.7071067811865475*uy[2]; 
  uz_xr[0] = 1.224744871391589*uz[1]+0.7071067811865475*uz[0]; 
  uz_xr[1] = 1.224744871391589*uz[3]+0.7071067811865475*uz[2]; 
 
  double *ux_yl = &u_surf[12]; 
  double *ux_yr = &u_surf[14]; 
  double *uy_yl = &u_surf[16]; 
  double *uy_yr = &u_surf[18]; 
  double *uz_yl = &u_surf[20]; 
  double *uz_yr = &u_surf[22]; 
 
  ux_yl[0] = 0.7071067811865475*ux[0]-1.224744871391589*ux[2]; 
  ux_yl[1] = 0.7071067811865475*ux[1]-1.224744871391589*ux[3]; 
  uy_yl[0] = 0.7071067811865475*uy[0]-1.224744871391589*uy[2]; 
  uy_yl[1] = 0.7071067811865475*uy[1]-1.224744871391589*uy[3]; 
  uz_yl[0] = 0.7071067811865475*uz[0]-1.224744871391589*uz[2]; 
  uz_yl[1] = 0.7071067811865475*uz[1]-1.224744871391589*uz[3]; 
 
  ux_yr[0] = 1.224744871391589*ux[2]+0.7071067811865475*ux[0]; 
  ux_yr[1] = 1.224744871391589*ux[3]+0.7071067811865475*ux[1]; 
  uy_yr[0] = 1.224744871391589*uy[2]+0.7071067811865475*uy[0]; 
  uy_yr[1] = 1.224744871391589*uy[3]+0.7071067811865475*uy[1]; 
  uz_yr[0] = 1.224744871391589*uz[2]+0.7071067811865475*uz[0]; 
  uz_yr[1] = 1.224744871391589*uz[3]+0.7071067811865475*uz[1]; 
 
} 
 
