#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH double lbo_gyrokinetic_drag_vol_1x1v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_gyrokinetic_diff_vol_1x1v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_gyrokinetic_drag_vol_1x1v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_boundary_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_gyrokinetic_diff_vol_1x1v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_gyrokinetic_drag_vol_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_boundary_surfvpar_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfvpar_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_boundary_surfmu_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfmu_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_gyrokinetic_diff_vol_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfvpar_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfvpar_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfmu_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfmu_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_gyrokinetic_drag_vol_1x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_boundary_surfvpar_1x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfvpar_1x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_boundary_surfmu_1x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfmu_1x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_gyrokinetic_diff_vol_1x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfvpar_1x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfvpar_1x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfmu_1x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfmu_1x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_gyrokinetic_drag_vol_2x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_boundary_surfvpar_2x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfvpar_2x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_boundary_surfmu_2x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfmu_2x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_gyrokinetic_diff_vol_2x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfvpar_2x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfvpar_2x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfmu_2x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfmu_2x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_gyrokinetic_drag_vol_2x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_boundary_surfvpar_2x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfvpar_2x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_boundary_surfmu_2x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfmu_2x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_gyrokinetic_diff_vol_2x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfvpar_2x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfvpar_2x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfmu_2x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfmu_2x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_gyrokinetic_drag_vol_3x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_boundary_surfvpar_3x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfvpar_3x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_boundary_surfmu_3x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfmu_3x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_gyrokinetic_diff_vol_3x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfvpar_3x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfvpar_3x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfmu_3x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfmu_3x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_gyrokinetic_drag_vol_3x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_boundary_surfvpar_3x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfvpar_3x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_boundary_surfmu_3x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfmu_3x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double lbo_gyrokinetic_diff_vol_3x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfvpar_3x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfvpar_3x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_boundary_surfmu_3x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out); 
GKYL_CU_DH double lbo_gyrokinetic_diff_surfmu_3x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 

EXTERN_C_END 
