#pragma once 
#include <math.h>
#include <gkyl_util.h>

EXTERN_C_BEG

typedef struct { double c, chi, gamma; } gkyl_maxwell_inp; 

GKYL_CU_DH double maxwell_vol_1x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_1x_ser_p1(const gkyl_maxwell_inp *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* GKYL_RESTRICT outl, double* GKYL_RESTRICT outr); 

GKYL_CU_DH double maxwell_vol_1x_ser_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_1x_ser_p2(const gkyl_maxwell_inp *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* GKYL_RESTRICT outl, double* GKYL_RESTRICT outr); 

GKYL_CU_DH double maxwell_vol_2x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_2x_ser_p1(const gkyl_maxwell_inp *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* GKYL_RESTRICT outl, double* GKYL_RESTRICT outr); 
GKYL_CU_DH double maxwell_surfy_2x_ser_p1(const gkyl_maxwell_inp *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* GKYL_RESTRICT outl, double* GKYL_RESTRICT outr); 

GKYL_CU_DH double maxwell_vol_2x_ser_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_2x_ser_p2(const gkyl_maxwell_inp *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* GKYL_RESTRICT outl, double* GKYL_RESTRICT outr); 
GKYL_CU_DH double maxwell_surfy_2x_ser_p2(const gkyl_maxwell_inp *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* GKYL_RESTRICT outl, double* GKYL_RESTRICT outr); 

GKYL_CU_DH double maxwell_vol_3x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_3x_ser_p1(const gkyl_maxwell_inp *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* GKYL_RESTRICT outl, double* GKYL_RESTRICT outr); 
GKYL_CU_DH double maxwell_surfy_3x_ser_p1(const gkyl_maxwell_inp *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* GKYL_RESTRICT outl, double* GKYL_RESTRICT outr); 
GKYL_CU_DH double maxwell_surfz_3x_ser_p1(const gkyl_maxwell_inp *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* GKYL_RESTRICT outl, double* GKYL_RESTRICT outr); 

GKYL_CU_DH double maxwell_vol_1x_tensor_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_1x_tensor_p2(const gkyl_maxwell_inp *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* GKYL_RESTRICT outl, double* GKYL_RESTRICT outr); 

GKYL_CU_DH double maxwell_vol_2x_tensor_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out); 
GKYL_CU_DH double maxwell_surfx_2x_tensor_p2(const gkyl_maxwell_inp *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* GKYL_RESTRICT outl, double* GKYL_RESTRICT outr); 
GKYL_CU_DH double maxwell_surfy_2x_tensor_p2(const gkyl_maxwell_inp *meq, const double *wl, const double *wr, const double *dxl, const double *dxr, const double tau, const double *ql, const double *qr, double* GKYL_RESTRICT outl, double* GKYL_RESTRICT outr); 

EXTERN_C_END
