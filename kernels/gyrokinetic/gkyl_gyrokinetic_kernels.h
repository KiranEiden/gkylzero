#pragma once

#include <math.h>
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH double gyrokinetic_vol_1x1v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double* dApardt, const double *f, double* GKYL_RESTRICT out); 
  double gyrokinetic_vol_step2_1x1v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfx_1x1v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfvpar_1x1v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x1v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const int edge, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_1x1v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double* dApardt, const double *f, double* GKYL_RESTRICT out); 
  double gyrokinetic_vol_step2_1x1v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfx_1x1v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfvpar_1x1v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x1v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const int edge, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_1x2v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double* dApardt, const double *f, double* GKYL_RESTRICT out); 
  double gyrokinetic_vol_step2_1x2v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfx_1x2v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfvpar_1x2v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x2v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const int edge, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_1x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double* dApardt, const double *f, double* GKYL_RESTRICT out); 
  double gyrokinetic_vol_step2_1x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfx_1x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfvpar_1x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const int edge, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_2x2v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double* dApardt, const double *f, double* GKYL_RESTRICT out); 
  double gyrokinetic_vol_step2_2x2v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfx_2x2v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfy_2x2v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfvpar_2x2v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_2x2v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const int edge, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_2x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double* dApardt, const double *f, double* GKYL_RESTRICT out); 
  double gyrokinetic_vol_step2_2x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfx_2x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfy_2x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfvpar_2x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_2x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const int edge, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_3x2v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double* dApardt, const double *f, double* GKYL_RESTRICT out); 
  double gyrokinetic_vol_step2_3x2v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfx_3x2v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfy_3x2v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfz_3x2v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfvpar_3x2v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_3x2v_ser_p1(const double q_, const double m_, const double *w, const double *dxv, const int edge, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_3x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double* dApardt, const double *f, double* GKYL_RESTRICT out); 
  double gyrokinetic_vol_step2_3x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfx_3x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfy_3x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfz_3x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfvpar_3x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_3x2v_ser_p2(const double q_, const double m_, const double *w, const double *dxv, const int edge, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_1x1v_tensor_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double* dApardt, const double *f, double* GKYL_RESTRICT out); 
  double gyrokinetic_vol_step2_1x1v_tensor_p1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfx_1x1v_tensor_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfvpar_1x1v_tensor_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x1v_tensor_p1(const double q_, const double m_, const double *w, const double *dxv, const int edge, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_1x1v_tensor_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double* dApardt, const double *f, double* GKYL_RESTRICT out); 
  double gyrokinetic_vol_step2_1x1v_tensor_p2(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfx_1x1v_tensor_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfvpar_1x1v_tensor_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x1v_tensor_p2(const double q_, const double m_, const double *w, const double *dxv, const int edge, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_1x2v_tensor_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double* dApardt, const double *f, double* GKYL_RESTRICT out); 
  double gyrokinetic_vol_step2_1x2v_tensor_p1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfx_1x2v_tensor_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfvpar_1x2v_tensor_p1(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x2v_tensor_p1(const double q_, const double m_, const double *w, const double *dxv, const int edge, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 

GKYL_CU_DH double gyrokinetic_vol_1x2v_tensor_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double* dApardt, const double *f, double* GKYL_RESTRICT out); 
  double gyrokinetic_vol_step2_1x2v_tensor_p2(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfx_1x2v_tensor_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_surfvpar_1x2v_tensor_p2(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x2v_tensor_p2(const double q_, const double m_, const double *w, const double *dxv, const int edge, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_i, const double *phi, const double *Apar, const double *dApardt, const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out); 


EXTERN_C_END
