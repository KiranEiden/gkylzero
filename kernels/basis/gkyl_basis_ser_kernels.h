// Wed Oct  6 12:06:43 2021
#pragma once
#include <gkyl_util.h>
EXTERN_C_BEG
GKYL_CU_DH void eval_1d_ser_p0(const double *z, double *b);
GKYL_CU_DH double eval_expand_1d_ser_p0(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_1d_ser_p0(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_1d_ser_p0(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_1d_ser_p0(double *node_coords);
GKYL_CU_DH void nodal_to_modal_1d_ser_p0(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_1d_ser_p1(const double *z, double *b);
GKYL_CU_DH double eval_expand_1d_ser_p1(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_1d_ser_p1(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_1d_ser_p1(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_1d_ser_p1(double *node_coords);
GKYL_CU_DH void nodal_to_modal_1d_ser_p1(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_1d_ser_p2(const double *z, double *b);
GKYL_CU_DH double eval_expand_1d_ser_p2(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_1d_ser_p2(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_1d_ser_p2(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_1d_ser_p2(double *node_coords);
GKYL_CU_DH void nodal_to_modal_1d_ser_p2(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_1d_ser_p3(const double *z, double *b);
GKYL_CU_DH double eval_expand_1d_ser_p3(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_1d_ser_p3(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_1d_ser_p3(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_1d_ser_p3(double *node_coords);
GKYL_CU_DH void nodal_to_modal_1d_ser_p3(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_2d_ser_p0(const double *z, double *b);
GKYL_CU_DH double eval_expand_2d_ser_p0(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_2d_ser_p0(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_2d_ser_p0(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_2d_ser_p0(double *node_coords);
GKYL_CU_DH void nodal_to_modal_2d_ser_p0(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_2d_ser_p1(const double *z, double *b);
GKYL_CU_DH double eval_expand_2d_ser_p1(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_2d_ser_p1(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_2d_ser_p1(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_2d_ser_p1(double *node_coords);
GKYL_CU_DH void nodal_to_modal_2d_ser_p1(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_2d_ser_p2(const double *z, double *b);
GKYL_CU_DH double eval_expand_2d_ser_p2(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_2d_ser_p2(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_2d_ser_p2(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_2d_ser_p2(double *node_coords);
GKYL_CU_DH void nodal_to_modal_2d_ser_p2(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_2d_ser_p3(const double *z, double *b);
GKYL_CU_DH double eval_expand_2d_ser_p3(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_2d_ser_p3(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_2d_ser_p3(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_2d_ser_p3(double *node_coords);
GKYL_CU_DH void nodal_to_modal_2d_ser_p3(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_3d_ser_p0(const double *z, double *b);
GKYL_CU_DH double eval_expand_3d_ser_p0(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_3d_ser_p0(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_3d_ser_p0(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_3d_ser_p0(double *node_coords);
GKYL_CU_DH void nodal_to_modal_3d_ser_p0(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_3d_ser_p1(const double *z, double *b);
GKYL_CU_DH double eval_expand_3d_ser_p1(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_3d_ser_p1(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_3d_ser_p1(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_3d_ser_p1(double *node_coords);
GKYL_CU_DH void nodal_to_modal_3d_ser_p1(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_3d_ser_p2(const double *z, double *b);
GKYL_CU_DH double eval_expand_3d_ser_p2(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_3d_ser_p2(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_3d_ser_p2(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_3d_ser_p2(double *node_coords);
GKYL_CU_DH void nodal_to_modal_3d_ser_p2(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_3d_ser_p3(const double *z, double *b);
GKYL_CU_DH double eval_expand_3d_ser_p3(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_3d_ser_p3(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_3d_ser_p3(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_3d_ser_p3(double *node_coords);
GKYL_CU_DH void nodal_to_modal_3d_ser_p3(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_4d_ser_p0(const double *z, double *b);
GKYL_CU_DH double eval_expand_4d_ser_p0(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_4d_ser_p0(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_4d_ser_p0(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_4d_ser_p0(double *node_coords);
GKYL_CU_DH void nodal_to_modal_4d_ser_p0(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_4d_ser_p1(const double *z, double *b);
GKYL_CU_DH double eval_expand_4d_ser_p1(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_4d_ser_p1(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_4d_ser_p1(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_4d_ser_p1(double *node_coords);
GKYL_CU_DH void nodal_to_modal_4d_ser_p1(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_4d_ser_p2(const double *z, double *b);
GKYL_CU_DH double eval_expand_4d_ser_p2(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_4d_ser_p2(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_4d_ser_p2(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_4d_ser_p2(double *node_coords);
GKYL_CU_DH void nodal_to_modal_4d_ser_p2(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_4d_ser_p3(const double *z, double *b);
GKYL_CU_DH double eval_expand_4d_ser_p3(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_4d_ser_p3(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_4d_ser_p3(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_4d_ser_p3(double *node_coords);
GKYL_CU_DH void nodal_to_modal_4d_ser_p3(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_5d_ser_p0(const double *z, double *b);
GKYL_CU_DH double eval_expand_5d_ser_p0(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_5d_ser_p0(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_5d_ser_p0(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_5d_ser_p0(double *node_coords);
GKYL_CU_DH void nodal_to_modal_5d_ser_p0(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_5d_ser_p1(const double *z, double *b);
GKYL_CU_DH double eval_expand_5d_ser_p1(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_5d_ser_p1(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_5d_ser_p1(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_5d_ser_p1(double *node_coords);
GKYL_CU_DH void nodal_to_modal_5d_ser_p1(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_5d_ser_p2(const double *z, double *b);
GKYL_CU_DH double eval_expand_5d_ser_p2(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_5d_ser_p2(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_5d_ser_p2(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_5d_ser_p2(double *node_coords);
GKYL_CU_DH void nodal_to_modal_5d_ser_p2(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_5d_ser_p3(const double *z, double *b);
GKYL_CU_DH double eval_expand_5d_ser_p3(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_5d_ser_p3(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_5d_ser_p3(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_5d_ser_p3(double *node_coords);
GKYL_CU_DH void nodal_to_modal_5d_ser_p3(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_6d_ser_p0(const double *z, double *b);
GKYL_CU_DH double eval_expand_6d_ser_p0(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_6d_ser_p0(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_6d_ser_p0(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_6d_ser_p0(double *node_coords);
GKYL_CU_DH void nodal_to_modal_6d_ser_p0(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_6d_ser_p1(const double *z, double *b);
GKYL_CU_DH double eval_expand_6d_ser_p1(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_6d_ser_p1(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_6d_ser_p1(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_6d_ser_p1(double *node_coords);
GKYL_CU_DH void nodal_to_modal_6d_ser_p1(const double *fnodal, double *fmodal);
GKYL_CU_DH void eval_6d_ser_p2(const double *z, double *b);
GKYL_CU_DH double eval_expand_6d_ser_p2(const double *z, const double *f);
GKYL_CU_DH double eval_grad_expand_6d_ser_p2(int dir, const double *z, const double *f);
GKYL_CU_DH void flip_sign_6d_ser_p2(int dir, const double *f, double *fout );
GKYL_CU_DH void node_coords_6d_ser_p2(double *node_coords);
GKYL_CU_DH void nodal_to_modal_6d_ser_p2(const double *fnodal, double *fmodal);
EXTERN_C_END
