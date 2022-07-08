#include <gkyl_dg_gen_diffusion_kernels.h>

GKYL_CU_DH void
dg_gen_diffusion_surfyy_3x_ser_p1(const double* w, const double* dx,
  const double* Dij, const double* q[], double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // Dij: Diffusion coefficient in the center cell
  // q: Input field in the left cell
  // out: Incremented output

  const double Jyy = 4/dx[1]/dx[1];
  const double* Dyy = &Dij[24];

  const double* qclc = q[10];
  const double* qccc = q[13];
  const double* qcuc = q[16];

  out[0] += Jyy*((-0.3314563036811939*Dyy[7]*qcuc[7])-0.1913663861549356*Dyy[5]*qcuc[7]-0.3314563036811939*Dyy[7]*qclc[7]+0.1913663861549356*Dyy[5]*qclc[7]-0.6629126073623879*Dyy[7]*qccc[7]+0.3444594950788841*qcuc[5]*Dyy[7]-0.3444594950788841*qclc[5]*Dyy[7]-0.3314563036811939*Dyy[6]*qcuc[6]-0.1913663861549356*Dyy[3]*qcuc[6]-0.3314563036811939*Dyy[6]*qclc[6]+0.1913663861549356*Dyy[3]*qclc[6]-0.6629126073623879*Dyy[6]*qccc[6]+0.3444594950788841*qcuc[3]*Dyy[6]-0.3444594950788841*qclc[3]*Dyy[6]+0.1988737822087163*Dyy[5]*qcuc[5]+0.1988737822087163*Dyy[5]*qclc[5]-0.3977475644174328*Dyy[5]*qccc[5]-0.3314563036811939*Dyy[4]*qcuc[4]-0.1913663861549356*Dyy[1]*qcuc[4]-0.3314563036811939*Dyy[4]*qclc[4]+0.1913663861549356*Dyy[1]*qclc[4]-0.6629126073623879*Dyy[4]*qccc[4]+0.3444594950788841*qcuc[1]*Dyy[4]-0.3444594950788841*qclc[1]*Dyy[4]+0.1988737822087163*Dyy[3]*qcuc[3]+0.1988737822087163*Dyy[3]*qclc[3]-0.3977475644174328*Dyy[3]*qccc[3]-0.3314563036811939*Dyy[2]*qcuc[2]-0.1913663861549356*Dyy[0]*qcuc[2]-0.3314563036811939*Dyy[2]*qclc[2]+0.1913663861549356*Dyy[0]*qclc[2]-0.6629126073623879*Dyy[2]*qccc[2]+0.3444594950788841*qcuc[0]*Dyy[2]-0.3444594950788841*qclc[0]*Dyy[2]+0.1988737822087163*Dyy[1]*qcuc[1]+0.1988737822087163*Dyy[1]*qclc[1]-0.3977475644174328*Dyy[1]*qccc[1]+0.1988737822087163*Dyy[0]*qcuc[0]+0.1988737822087163*Dyy[0]*qclc[0]-0.3977475644174328*Dyy[0]*qccc[0]);
  out[1] += Jyy*((-0.3314563036811939*Dyy[6]*qcuc[7])-0.1913663861549356*Dyy[3]*qcuc[7]-0.3314563036811939*Dyy[6]*qclc[7]+0.1913663861549356*Dyy[3]*qclc[7]-0.6629126073623879*Dyy[6]*qccc[7]-0.3314563036811939*qcuc[6]*Dyy[7]-0.3314563036811939*qclc[6]*Dyy[7]-0.6629126073623879*qccc[6]*Dyy[7]+0.3444594950788841*qcuc[3]*Dyy[7]-0.3444594950788841*qclc[3]*Dyy[7]-0.1913663861549356*Dyy[5]*qcuc[6]+0.1913663861549356*Dyy[5]*qclc[6]+0.3444594950788841*qcuc[5]*Dyy[6]-0.3444594950788841*qclc[5]*Dyy[6]+0.1988737822087163*Dyy[3]*qcuc[5]+0.1988737822087163*Dyy[3]*qclc[5]-0.3977475644174328*Dyy[3]*qccc[5]+0.1988737822087163*qcuc[3]*Dyy[5]+0.1988737822087163*qclc[3]*Dyy[5]-0.3977475644174328*qccc[3]*Dyy[5]-0.3314563036811939*Dyy[2]*qcuc[4]-0.1913663861549356*Dyy[0]*qcuc[4]-0.3314563036811939*Dyy[2]*qclc[4]+0.1913663861549356*Dyy[0]*qclc[4]-0.6629126073623879*Dyy[2]*qccc[4]-0.3314563036811939*qcuc[2]*Dyy[4]-0.3314563036811939*qclc[2]*Dyy[4]-0.6629126073623879*qccc[2]*Dyy[4]+0.3444594950788841*qcuc[0]*Dyy[4]-0.3444594950788841*qclc[0]*Dyy[4]-0.1913663861549356*Dyy[1]*qcuc[2]+0.1913663861549356*Dyy[1]*qclc[2]+0.3444594950788841*qcuc[1]*Dyy[2]-0.3444594950788841*qclc[1]*Dyy[2]+0.1988737822087163*Dyy[0]*qcuc[1]+0.1988737822087163*Dyy[0]*qclc[1]-0.3977475644174328*Dyy[0]*qccc[1]+0.1988737822087163*qcuc[0]*Dyy[1]+0.1988737822087163*qclc[0]*Dyy[1]-0.3977475644174328*qccc[0]*Dyy[1]);
  out[2] += Jyy*((-0.2679129406169099*Dyy[7]*qcuc[7])-0.1546796083845572*Dyy[5]*qcuc[7]+0.2679129406169099*Dyy[7]*qclc[7]-0.1546796083845572*Dyy[5]*qclc[7]-1.016465997955661*Dyy[5]*qccc[7]+0.3314563036811939*qcuc[5]*Dyy[7]+0.3314563036811939*qclc[5]*Dyy[7]-0.6629126073623879*qccc[5]*Dyy[7]-0.2679129406169099*Dyy[6]*qcuc[6]-0.1546796083845572*Dyy[3]*qcuc[6]+0.2679129406169099*Dyy[6]*qclc[6]-0.1546796083845572*Dyy[3]*qclc[6]-1.016465997955661*Dyy[3]*qccc[6]+0.3314563036811939*qcuc[3]*Dyy[6]+0.3314563036811939*qclc[3]*Dyy[6]-0.6629126073623879*qccc[3]*Dyy[6]+0.1913663861549356*Dyy[5]*qcuc[5]-0.1913663861549356*Dyy[5]*qclc[5]-0.2679129406169099*Dyy[4]*qcuc[4]-0.1546796083845572*Dyy[1]*qcuc[4]+0.2679129406169099*Dyy[4]*qclc[4]-0.1546796083845572*Dyy[1]*qclc[4]-1.016465997955661*Dyy[1]*qccc[4]+0.3314563036811939*qcuc[1]*Dyy[4]+0.3314563036811939*qclc[1]*Dyy[4]-0.6629126073623879*qccc[1]*Dyy[4]+0.1913663861549356*Dyy[3]*qcuc[3]-0.1913663861549356*Dyy[3]*qclc[3]-0.2679129406169099*Dyy[2]*qcuc[2]-0.1546796083845572*Dyy[0]*qcuc[2]+0.2679129406169099*Dyy[2]*qclc[2]-0.1546796083845572*Dyy[0]*qclc[2]-1.016465997955661*Dyy[0]*qccc[2]+0.3314563036811939*qcuc[0]*Dyy[2]+0.3314563036811939*qclc[0]*Dyy[2]-0.6629126073623879*qccc[0]*Dyy[2]+0.1913663861549356*Dyy[1]*qcuc[1]-0.1913663861549356*Dyy[1]*qclc[1]+0.1913663861549356*Dyy[0]*qcuc[0]-0.1913663861549356*Dyy[0]*qclc[0]);
  out[3] += Jyy*((-0.3314563036811939*Dyy[4]*qcuc[7])-0.1913663861549356*Dyy[1]*qcuc[7]-0.3314563036811939*Dyy[4]*qclc[7]+0.1913663861549356*Dyy[1]*qclc[7]-0.6629126073623879*Dyy[4]*qccc[7]-0.3314563036811939*qcuc[4]*Dyy[7]-0.3314563036811939*qclc[4]*Dyy[7]-0.6629126073623879*qccc[4]*Dyy[7]+0.3444594950788841*qcuc[1]*Dyy[7]-0.3444594950788841*qclc[1]*Dyy[7]-0.3314563036811939*Dyy[2]*qcuc[6]-0.1913663861549356*Dyy[0]*qcuc[6]-0.3314563036811939*Dyy[2]*qclc[6]+0.1913663861549356*Dyy[0]*qclc[6]-0.6629126073623879*Dyy[2]*qccc[6]-0.3314563036811939*qcuc[2]*Dyy[6]-0.3314563036811939*qclc[2]*Dyy[6]-0.6629126073623879*qccc[2]*Dyy[6]+0.3444594950788841*qcuc[0]*Dyy[6]-0.3444594950788841*qclc[0]*Dyy[6]+0.3444594950788841*Dyy[4]*qcuc[5]+0.1988737822087163*Dyy[1]*qcuc[5]-0.3444594950788841*Dyy[4]*qclc[5]+0.1988737822087163*Dyy[1]*qclc[5]-0.3977475644174328*Dyy[1]*qccc[5]-0.1913663861549356*qcuc[4]*Dyy[5]+0.1913663861549356*qclc[4]*Dyy[5]+0.1988737822087163*qcuc[1]*Dyy[5]+0.1988737822087163*qclc[1]*Dyy[5]-0.3977475644174328*qccc[1]*Dyy[5]+0.3444594950788841*Dyy[2]*qcuc[3]+0.1988737822087163*Dyy[0]*qcuc[3]-0.3444594950788841*Dyy[2]*qclc[3]+0.1988737822087163*Dyy[0]*qclc[3]-0.3977475644174328*Dyy[0]*qccc[3]-0.1913663861549356*qcuc[2]*Dyy[3]+0.1913663861549356*qclc[2]*Dyy[3]+0.1988737822087163*qcuc[0]*Dyy[3]+0.1988737822087163*qclc[0]*Dyy[3]-0.3977475644174328*qccc[0]*Dyy[3]);
  out[4] += Jyy*((-0.2679129406169099*Dyy[6]*qcuc[7])-0.1546796083845572*Dyy[3]*qcuc[7]+0.2679129406169099*Dyy[6]*qclc[7]-0.1546796083845572*Dyy[3]*qclc[7]-1.016465997955661*Dyy[3]*qccc[7]-0.2679129406169099*qcuc[6]*Dyy[7]+0.2679129406169099*qclc[6]*Dyy[7]+0.3314563036811939*qcuc[3]*Dyy[7]+0.3314563036811939*qclc[3]*Dyy[7]-0.6629126073623879*qccc[3]*Dyy[7]-0.1546796083845572*Dyy[5]*qcuc[6]-0.1546796083845572*Dyy[5]*qclc[6]-1.016465997955661*Dyy[5]*qccc[6]+0.3314563036811939*qcuc[5]*Dyy[6]+0.3314563036811939*qclc[5]*Dyy[6]-0.6629126073623879*qccc[5]*Dyy[6]+0.1913663861549356*Dyy[3]*qcuc[5]-0.1913663861549356*Dyy[3]*qclc[5]+0.1913663861549356*qcuc[3]*Dyy[5]-0.1913663861549356*qclc[3]*Dyy[5]-0.2679129406169099*Dyy[2]*qcuc[4]-0.1546796083845572*Dyy[0]*qcuc[4]+0.2679129406169099*Dyy[2]*qclc[4]-0.1546796083845572*Dyy[0]*qclc[4]-1.016465997955661*Dyy[0]*qccc[4]-0.2679129406169099*qcuc[2]*Dyy[4]+0.2679129406169099*qclc[2]*Dyy[4]+0.3314563036811939*qcuc[0]*Dyy[4]+0.3314563036811939*qclc[0]*Dyy[4]-0.6629126073623879*qccc[0]*Dyy[4]-0.1546796083845572*Dyy[1]*qcuc[2]-0.1546796083845572*Dyy[1]*qclc[2]-1.016465997955661*Dyy[1]*qccc[2]+0.3314563036811939*qcuc[1]*Dyy[2]+0.3314563036811939*qclc[1]*Dyy[2]-0.6629126073623879*qccc[1]*Dyy[2]+0.1913663861549356*Dyy[0]*qcuc[1]-0.1913663861549356*Dyy[0]*qclc[1]+0.1913663861549356*qcuc[0]*Dyy[1]-0.1913663861549356*qclc[0]*Dyy[1]);
  out[5] += Jyy*((-0.3314563036811939*Dyy[2]*qcuc[7])-0.1913663861549356*Dyy[0]*qcuc[7]-0.3314563036811939*Dyy[2]*qclc[7]+0.1913663861549356*Dyy[0]*qclc[7]-0.6629126073623879*Dyy[2]*qccc[7]-0.3314563036811939*qcuc[2]*Dyy[7]-0.3314563036811939*qclc[2]*Dyy[7]-0.6629126073623879*qccc[2]*Dyy[7]+0.3444594950788841*qcuc[0]*Dyy[7]-0.3444594950788841*qclc[0]*Dyy[7]-0.3314563036811939*Dyy[4]*qcuc[6]-0.1913663861549356*Dyy[1]*qcuc[6]-0.3314563036811939*Dyy[4]*qclc[6]+0.1913663861549356*Dyy[1]*qclc[6]-0.6629126073623879*Dyy[4]*qccc[6]-0.3314563036811939*qcuc[4]*Dyy[6]-0.3314563036811939*qclc[4]*Dyy[6]-0.6629126073623879*qccc[4]*Dyy[6]+0.3444594950788841*qcuc[1]*Dyy[6]-0.3444594950788841*qclc[1]*Dyy[6]+0.3444594950788841*Dyy[2]*qcuc[5]+0.1988737822087163*Dyy[0]*qcuc[5]-0.3444594950788841*Dyy[2]*qclc[5]+0.1988737822087163*Dyy[0]*qclc[5]-0.3977475644174328*Dyy[0]*qccc[5]-0.1913663861549356*qcuc[2]*Dyy[5]+0.1913663861549356*qclc[2]*Dyy[5]+0.1988737822087163*qcuc[0]*Dyy[5]+0.1988737822087163*qclc[0]*Dyy[5]-0.3977475644174328*qccc[0]*Dyy[5]-0.1913663861549356*Dyy[3]*qcuc[4]+0.1913663861549356*Dyy[3]*qclc[4]+0.3444594950788841*qcuc[3]*Dyy[4]-0.3444594950788841*qclc[3]*Dyy[4]+0.1988737822087163*Dyy[1]*qcuc[3]+0.1988737822087163*Dyy[1]*qclc[3]-0.3977475644174328*Dyy[1]*qccc[3]+0.1988737822087163*qcuc[1]*Dyy[3]+0.1988737822087163*qclc[1]*Dyy[3]-0.3977475644174328*qccc[1]*Dyy[3]);
  out[6] += Jyy*((-0.2679129406169099*Dyy[4]*qcuc[7])-0.1546796083845572*Dyy[1]*qcuc[7]+0.2679129406169099*Dyy[4]*qclc[7]-0.1546796083845572*Dyy[1]*qclc[7]-1.016465997955661*Dyy[1]*qccc[7]-0.2679129406169099*qcuc[4]*Dyy[7]+0.2679129406169099*qclc[4]*Dyy[7]+0.3314563036811939*qcuc[1]*Dyy[7]+0.3314563036811939*qclc[1]*Dyy[7]-0.6629126073623879*qccc[1]*Dyy[7]-0.2679129406169099*Dyy[2]*qcuc[6]-0.1546796083845572*Dyy[0]*qcuc[6]+0.2679129406169099*Dyy[2]*qclc[6]-0.1546796083845572*Dyy[0]*qclc[6]-1.016465997955661*Dyy[0]*qccc[6]-0.2679129406169099*qcuc[2]*Dyy[6]+0.2679129406169099*qclc[2]*Dyy[6]+0.3314563036811939*qcuc[0]*Dyy[6]+0.3314563036811939*qclc[0]*Dyy[6]-0.6629126073623879*qccc[0]*Dyy[6]+0.3314563036811939*Dyy[4]*qcuc[5]+0.1913663861549356*Dyy[1]*qcuc[5]+0.3314563036811939*Dyy[4]*qclc[5]-0.1913663861549356*Dyy[1]*qclc[5]-0.6629126073623879*Dyy[4]*qccc[5]-0.1546796083845572*qcuc[4]*Dyy[5]-0.1546796083845572*qclc[4]*Dyy[5]-1.016465997955661*qccc[4]*Dyy[5]+0.1913663861549356*qcuc[1]*Dyy[5]-0.1913663861549356*qclc[1]*Dyy[5]+0.3314563036811939*Dyy[2]*qcuc[3]+0.1913663861549356*Dyy[0]*qcuc[3]+0.3314563036811939*Dyy[2]*qclc[3]-0.1913663861549356*Dyy[0]*qclc[3]-0.6629126073623879*Dyy[2]*qccc[3]-0.1546796083845572*qcuc[2]*Dyy[3]-0.1546796083845572*qclc[2]*Dyy[3]-1.016465997955661*qccc[2]*Dyy[3]+0.1913663861549356*qcuc[0]*Dyy[3]-0.1913663861549356*qclc[0]*Dyy[3]);
  out[7] += Jyy*((-0.2679129406169099*Dyy[2]*qcuc[7])-0.1546796083845572*Dyy[0]*qcuc[7]+0.2679129406169099*Dyy[2]*qclc[7]-0.1546796083845572*Dyy[0]*qclc[7]-1.016465997955661*Dyy[0]*qccc[7]-0.2679129406169099*qcuc[2]*Dyy[7]+0.2679129406169099*qclc[2]*Dyy[7]+0.3314563036811939*qcuc[0]*Dyy[7]+0.3314563036811939*qclc[0]*Dyy[7]-0.6629126073623879*qccc[0]*Dyy[7]-0.2679129406169099*Dyy[4]*qcuc[6]-0.1546796083845572*Dyy[1]*qcuc[6]+0.2679129406169099*Dyy[4]*qclc[6]-0.1546796083845572*Dyy[1]*qclc[6]-1.016465997955661*Dyy[1]*qccc[6]-0.2679129406169099*qcuc[4]*Dyy[6]+0.2679129406169099*qclc[4]*Dyy[6]+0.3314563036811939*qcuc[1]*Dyy[6]+0.3314563036811939*qclc[1]*Dyy[6]-0.6629126073623879*qccc[1]*Dyy[6]+0.3314563036811939*Dyy[2]*qcuc[5]+0.1913663861549356*Dyy[0]*qcuc[5]+0.3314563036811939*Dyy[2]*qclc[5]-0.1913663861549356*Dyy[0]*qclc[5]-0.6629126073623879*Dyy[2]*qccc[5]-0.1546796083845572*qcuc[2]*Dyy[5]-0.1546796083845572*qclc[2]*Dyy[5]-1.016465997955661*qccc[2]*Dyy[5]+0.1913663861549356*qcuc[0]*Dyy[5]-0.1913663861549356*qclc[0]*Dyy[5]-0.1546796083845572*Dyy[3]*qcuc[4]-0.1546796083845572*Dyy[3]*qclc[4]-1.016465997955661*Dyy[3]*qccc[4]+0.3314563036811939*qcuc[3]*Dyy[4]+0.3314563036811939*qclc[3]*Dyy[4]-0.6629126073623879*qccc[3]*Dyy[4]+0.1913663861549356*Dyy[1]*qcuc[3]-0.1913663861549356*Dyy[1]*qclc[3]+0.1913663861549356*qcuc[1]*Dyy[3]-0.1913663861549356*qclc[1]*Dyy[3]);
}
