#pragma once

#include <math.h>
#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH double canonical_pb_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *hamil,  const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double *hamil, 
            double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double *hamil,
              double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double *hamil, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double *hamil,
          const double *alpha_surf_edge, const double *alpha_surf_skin, 
          const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
          const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
          const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *hamil, 
            double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *hamil, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *hamil,
          const double *alpha_surf_edge, const double *alpha_surf_skin, 
          const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
          const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
          const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *hamil,  const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfx_1x1v_ser_p2(const double *w, const double *dxv, const double *hamil, 
            double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_1x1v_ser_p2(const double *w, const double *dxv, const double *hamil,
              double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_1x1v_ser_p2(const double *w, const double *dxv, const double *hamil, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv, const double *hamil,
          const double *alpha_surf_edge, const double *alpha_surf_skin, 
          const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
          const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
          const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *hamil, 
            double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *hamil, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *hamil,
          const double *alpha_surf_edge, const double *alpha_surf_skin, 
          const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
          const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
          const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil,  const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
            double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil,
              double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil,
          const double *alpha_surf_edge, const double *alpha_surf_skin, 
          const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
          const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
          const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
            double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil,
              double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil,
          const double *alpha_surf_edge, const double *alpha_surf_skin, 
          const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
          const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
          const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
            double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil,
              double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil,
          const double *alpha_surf_edge, const double *alpha_surf_skin, 
          const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
          const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
          const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
            double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *hamil,
          const double *alpha_surf_edge, const double *alpha_surf_skin, 
          const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
          const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
          const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil,  const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfx_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
            double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil,
              double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil,
          const double *alpha_surf_edge, const double *alpha_surf_skin, 
          const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
          const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
          const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfy_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
            double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil,
              double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfy_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfy_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil,
          const double *alpha_surf_edge, const double *alpha_surf_skin, 
          const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
          const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
          const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
            double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil,
              double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil,
          const double *alpha_surf_edge, const double *alpha_surf_skin, 
          const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
          const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
          const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
            double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *hamil,
          const double *alpha_surf_edge, const double *alpha_surf_skin, 
          const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
          const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
          const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 

GKYL_CU_DH double canonical_pb_vol_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil,  const double *fin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfx_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
            double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfx_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
              double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfx_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfx_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
          const double *alpha_surf_edge, const double *alpha_surf_skin, 
          const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
          const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
          const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfy_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
            double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfy_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
              double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfy_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfy_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
          const double *alpha_surf_edge, const double *alpha_surf_skin, 
          const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
          const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
          const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfz_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
            double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfz_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
              double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfz_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfz_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
          const double *alpha_surf_edge, const double *alpha_surf_skin, 
          const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
          const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
          const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
            double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
              double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
          const double *alpha_surf_edge, const double *alpha_surf_skin, 
          const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
          const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
          const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
            double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH int canonical_pb_alpha_edge_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
              double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
          const double *alpha_surf_edge, const double *alpha_surf_skin, 
          const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
          const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
          const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 
GKYL_CU_DH int canonical_pb_alpha_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
            double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
GKYL_CU_DH double canonical_pb_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil, 
          const double *alpha_surf_l, const double *alpha_surf_r, 
          const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
          const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
          const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out); 
GKYL_CU_DH double canonical_pb_boundary_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double *hamil,
          const double *alpha_surf_edge, const double *alpha_surf_skin, 
          const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
          const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
          const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out); 


EXTERN_C_END
