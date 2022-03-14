#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
#include <gkyl_util.h> 
 
EXTERN_C_BEG 
 
GKYL_CU_DH void prim_lbo_copy_sol(const struct gkyl_mat *rhs, const int nc, const int vdim, double* GKYL_RESTRICT u, double* GKYL_RESTRICT vtSq); 
 
GKYL_CU_DH void vlasov_self_prim_moments_1x1v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x1v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double m_self,
  const double *u_self, const double *vtsq_self, const double m_other, const double *u_other,
  const double *vtsq_other, const double *nu, const double *moms, const double *boundary_corrections); 

GKYL_CU_DH void vlasov_self_prim_moments_1x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double m_self,
  const double *u_self, const double *vtsq_self, const double m_other, const double *u_other,
  const double *vtsq_other, const double *nu, const double *moms, const double *boundary_corrections); 

GKYL_CU_DH void vlasov_self_prim_moments_1x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double m_self,
  const double *u_self, const double *vtsq_self, const double m_other, const double *u_other,
  const double *vtsq_other, const double *nu, const double *moms, const double *boundary_corrections); 

GKYL_CU_DH void vlasov_self_prim_moments_2x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double m_self,
  const double *u_self, const double *vtsq_self, const double m_other, const double *u_other,
  const double *vtsq_other, const double *nu, const double *moms, const double *boundary_corrections); 

GKYL_CU_DH void vlasov_self_prim_moments_2x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double m_self,
  const double *u_self, const double *vtsq_self, const double m_other, const double *u_other,
  const double *vtsq_other, const double *nu, const double *moms, const double *boundary_corrections); 

GKYL_CU_DH void vlasov_self_prim_moments_3x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *moms, const double *boundary_corrections); 
GKYL_CU_DH void vlasov_cross_prim_moments_3x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double m_self,
  const double *u_self, const double *vtsq_self, const double m_other, const double *u_other,
  const double *vtsq_other, const double *nu, const double *moms, const double *boundary_corrections); 


EXTERN_C_END 
