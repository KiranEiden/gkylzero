#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH void em_surf_set_bvar_1x_ser_p1(const double* bvar, double* GKYL_RESTRICT bvar_surf) 
{ 
  // bvar:      Input volume expansion of b_i = B_i/|B| (first 3 components), b_i b_j = B_i B_j/|B|^2 (last 6 components). 
  // bvar_surf: Output surface expansion of magnetic field unit vector and tensor in each direction. 
 
  const double *bx = &bvar[0]; 
  const double *by = &bvar[2]; 
  const double *bz = &bvar[4]; 
  const double *bxbx = &bvar[6]; 
  const double *bxby = &bvar[8]; 
  const double *bxbz = &bvar[10]; 
  const double *byby = &bvar[12]; 
  const double *bybz = &bvar[14]; 
  const double *bzbz = &bvar[16]; 
 
  double *bx_xl = &bvar_surf[0]; 
  double *bx_xr = &bvar_surf[1]; 
  double *bxbx_xl = &bvar_surf[2]; 
  double *bxbx_xr = &bvar_surf[3]; 
  double *bxby_xl = &bvar_surf[4]; 
  double *bxby_xr = &bvar_surf[5]; 
  double *bxbz_xl = &bvar_surf[6]; 
  double *bxbz_xr = &bvar_surf[7]; 
 
  bx_xl[0] = 0.7071067811865475*bx[0]-1.224744871391589*bx[1]; 
  bx_xr[0] = 1.224744871391589*bx[1]+0.7071067811865475*bx[0]; 
  bxbx_xl[0] = 0.7071067811865475*bxbx[0]-1.224744871391589*bxbx[1]; 
  bxbx_xr[0] = 1.224744871391589*bxbx[1]+0.7071067811865475*bxbx[0]; 
  bxby_xl[0] = 0.7071067811865475*bxby[0]-1.224744871391589*bxby[1]; 
  bxby_xr[0] = 1.224744871391589*bxby[1]+0.7071067811865475*bxby[0]; 
  bxbz_xl[0] = 0.7071067811865475*bxbz[0]-1.224744871391589*bxbz[1]; 
  bxbz_xr[0] = 1.224744871391589*bxbz[1]+0.7071067811865475*bxbz[0]; 
 
} 
 
