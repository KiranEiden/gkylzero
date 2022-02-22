#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void gyrokinetic_M0_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M1_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_par_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M3_par_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_ThreeMoments_1x1v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT outM0, double* GKYL_RESTRICT outM1, double* GKYL_RESTRICT outM2); 

GKYL_CU_DH void gyrokinetic_M0_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M1_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_par_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M3_par_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_ThreeMoments_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT outM0, double* GKYL_RESTRICT outM1, double* GKYL_RESTRICT outM2); 

GKYL_CU_DH void gyrokinetic_M0_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M1_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_par_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_perp_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M3_par_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M3_perp_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_ThreeMoments_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT outM0, double* GKYL_RESTRICT outM1, double* GKYL_RESTRICT outM2); 
GKYL_CU_DH void gyrokinetic_M0_step1_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M0_step2_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void gyrokinetic_M0_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M1_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_par_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_perp_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M3_par_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M3_perp_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_ThreeMoments_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT outM0, double* GKYL_RESTRICT outM1, double* GKYL_RESTRICT outM2); 
GKYL_CU_DH void gyrokinetic_M0_step1_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M0_step2_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void gyrokinetic_M0_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M1_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_par_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_perp_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M3_par_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M3_perp_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_ThreeMoments_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT outM0, double* GKYL_RESTRICT outM1, double* GKYL_RESTRICT outM2); 
GKYL_CU_DH void gyrokinetic_M0_step1_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M0_step2_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void gyrokinetic_M0_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M1_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_par_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_perp_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M3_par_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M3_perp_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_ThreeMoments_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT outM0, double* GKYL_RESTRICT outM1, double* GKYL_RESTRICT outM2); 
GKYL_CU_DH void gyrokinetic_M0_step1_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M0_step2_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void gyrokinetic_M0_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M1_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_par_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_perp_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M3_par_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M3_perp_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_ThreeMoments_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT outM0, double* GKYL_RESTRICT outM1, double* GKYL_RESTRICT outM2); 
GKYL_CU_DH void gyrokinetic_M0_step1_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M0_step2_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void gyrokinetic_M0_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M1_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_par_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M3_par_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_ThreeMoments_1x1v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT outM0, double* GKYL_RESTRICT outM1, double* GKYL_RESTRICT outM2); 

GKYL_CU_DH void gyrokinetic_M0_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M1_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_par_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_perp_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M3_par_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M3_perp_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_ThreeMoments_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT outM0, double* GKYL_RESTRICT outM1, double* GKYL_RESTRICT outM2); 
GKYL_CU_DH void gyrokinetic_M0_step1_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M0_step2_1x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 

GKYL_CU_DH void gyrokinetic_M0_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M1_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_par_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M2_perp_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M3_par_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M3_perp_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_ThreeMoments_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT outM0, double* GKYL_RESTRICT outM1, double* GKYL_RESTRICT outM2); 
GKYL_CU_DH void gyrokinetic_M0_step1_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH void gyrokinetic_M0_step2_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out); 

EXTERN_C_END 
