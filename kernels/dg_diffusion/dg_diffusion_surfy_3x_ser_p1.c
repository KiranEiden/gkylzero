#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH void
dg_diffusion_surfy_3x_ser_p1(const double* w, const double* dx,
  const double* D,
  const double* ql, const double* qc, const double* qr,
  double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // D: Diffusion coefficient in the center cell
  // ql: Input field in the left cell
  // qc: Input field in the center cell
  // qr: Input field in the right cell
  // out: Incremented output

  const double J = 4/dx[1]/dx[1];

  out[0] += J*((-0.3314563036811939*qr[7]*D[15])-0.3314563036811939*ql[7]*D[15]-0.6629126073623879*qc[7]*D[15]+0.3444594950788841*qr[5]*D[15]-0.3444594950788841*ql[5]*D[15]-0.3314563036811939*qr[6]*D[14]-0.3314563036811939*ql[6]*D[14]-0.6629126073623879*qc[6]*D[14]+0.3444594950788841*qr[3]*D[14]-0.3444594950788841*ql[3]*D[14]-0.1913663861549356*qr[7]*D[13]+0.1913663861549356*ql[7]*D[13]+0.1988737822087163*qr[5]*D[13]+0.1988737822087163*ql[5]*D[13]-0.3977475644174328*qc[5]*D[13]-0.3314563036811939*qr[4]*D[12]-0.3314563036811939*ql[4]*D[12]-0.6629126073623879*qc[4]*D[12]+0.3444594950788841*qr[1]*D[12]-0.3444594950788841*ql[1]*D[12]-0.1913663861549356*qr[6]*D[11]+0.1913663861549356*ql[6]*D[11]+0.1988737822087163*qr[3]*D[11]+0.1988737822087163*ql[3]*D[11]-0.3977475644174328*qc[3]*D[11]-0.3314563036811939*qr[2]*D[10]-0.3314563036811939*ql[2]*D[10]-0.6629126073623879*qc[2]*D[10]+0.3444594950788841*qr[0]*D[10]-0.3444594950788841*ql[0]*D[10]-0.1913663861549356*qr[4]*D[9]+0.1913663861549356*ql[4]*D[9]+0.1988737822087163*qr[1]*D[9]+0.1988737822087163*ql[1]*D[9]-0.3977475644174328*qc[1]*D[9]-0.1913663861549356*qr[2]*D[8]+0.1913663861549356*ql[2]*D[8]+0.1988737822087163*qr[0]*D[8]+0.1988737822087163*ql[0]*D[8]-0.3977475644174328*qc[0]*D[8]);
  out[1] += J*((-0.3314563036811939*qr[6]*D[15])-0.3314563036811939*ql[6]*D[15]-0.6629126073623879*qc[6]*D[15]+0.3444594950788841*qr[3]*D[15]-0.3444594950788841*ql[3]*D[15]-0.3314563036811939*qr[7]*D[14]-0.3314563036811939*ql[7]*D[14]-0.6629126073623879*qc[7]*D[14]+0.3444594950788841*qr[5]*D[14]-0.3444594950788841*ql[5]*D[14]-0.1913663861549356*qr[6]*D[13]+0.1913663861549356*ql[6]*D[13]+0.1988737822087163*qr[3]*D[13]+0.1988737822087163*ql[3]*D[13]-0.3977475644174328*qc[3]*D[13]-0.3314563036811939*qr[2]*D[12]-0.3314563036811939*ql[2]*D[12]-0.6629126073623879*qc[2]*D[12]+0.3444594950788841*qr[0]*D[12]-0.3444594950788841*ql[0]*D[12]-0.1913663861549356*qr[7]*D[11]+0.1913663861549356*ql[7]*D[11]+0.1988737822087163*qr[5]*D[11]+0.1988737822087163*ql[5]*D[11]-0.3977475644174328*qc[5]*D[11]-0.3314563036811939*qr[4]*D[10]-0.3314563036811939*ql[4]*D[10]-0.6629126073623879*qc[4]*D[10]+0.3444594950788841*qr[1]*D[10]-0.3444594950788841*ql[1]*D[10]-0.1913663861549356*qr[2]*D[9]+0.1913663861549356*ql[2]*D[9]+0.1988737822087163*qr[0]*D[9]+0.1988737822087163*ql[0]*D[9]-0.3977475644174328*qc[0]*D[9]-0.1913663861549356*qr[4]*D[8]+0.1913663861549356*ql[4]*D[8]+0.1988737822087163*qr[1]*D[8]+0.1988737822087163*ql[1]*D[8]-0.3977475644174328*qc[1]*D[8]);
  out[2] += J*((-0.2679129406169099*qr[7]*D[15])+0.2679129406169099*ql[7]*D[15]+0.3314563036811939*qr[5]*D[15]+0.3314563036811939*ql[5]*D[15]-0.6629126073623879*qc[5]*D[15]-0.2679129406169099*qr[6]*D[14]+0.2679129406169099*ql[6]*D[14]+0.3314563036811939*qr[3]*D[14]+0.3314563036811939*ql[3]*D[14]-0.6629126073623879*qc[3]*D[14]-0.1546796083845572*qr[7]*D[13]-0.1546796083845572*ql[7]*D[13]-1.016465997955661*qc[7]*D[13]+0.1913663861549356*qr[5]*D[13]-0.1913663861549356*ql[5]*D[13]-0.2679129406169099*qr[4]*D[12]+0.2679129406169099*ql[4]*D[12]+0.3314563036811939*qr[1]*D[12]+0.3314563036811939*ql[1]*D[12]-0.6629126073623879*qc[1]*D[12]-0.1546796083845572*qr[6]*D[11]-0.1546796083845572*ql[6]*D[11]-1.016465997955661*qc[6]*D[11]+0.1913663861549356*qr[3]*D[11]-0.1913663861549356*ql[3]*D[11]-0.2679129406169099*qr[2]*D[10]+0.2679129406169099*ql[2]*D[10]+0.3314563036811939*qr[0]*D[10]+0.3314563036811939*ql[0]*D[10]-0.6629126073623879*qc[0]*D[10]-0.1546796083845572*qr[4]*D[9]-0.1546796083845572*ql[4]*D[9]-1.016465997955661*qc[4]*D[9]+0.1913663861549356*qr[1]*D[9]-0.1913663861549356*ql[1]*D[9]-0.1546796083845572*qr[2]*D[8]-0.1546796083845572*ql[2]*D[8]-1.016465997955661*qc[2]*D[8]+0.1913663861549356*qr[0]*D[8]-0.1913663861549356*ql[0]*D[8]);
  out[3] += J*((-0.3314563036811939*qr[4]*D[15])-0.3314563036811939*ql[4]*D[15]-0.6629126073623879*qc[4]*D[15]+0.3444594950788841*qr[1]*D[15]-0.3444594950788841*ql[1]*D[15]-0.3314563036811939*qr[2]*D[14]-0.3314563036811939*ql[2]*D[14]-0.6629126073623879*qc[2]*D[14]+0.3444594950788841*qr[0]*D[14]-0.3444594950788841*ql[0]*D[14]-0.1913663861549356*qr[4]*D[13]+0.1913663861549356*ql[4]*D[13]+0.1988737822087163*qr[1]*D[13]+0.1988737822087163*ql[1]*D[13]-0.3977475644174328*qc[1]*D[13]-0.3314563036811939*qr[7]*D[12]-0.3314563036811939*ql[7]*D[12]-0.6629126073623879*qc[7]*D[12]+0.3444594950788841*qr[5]*D[12]-0.3444594950788841*ql[5]*D[12]-0.1913663861549356*qr[2]*D[11]+0.1913663861549356*ql[2]*D[11]+0.1988737822087163*qr[0]*D[11]+0.1988737822087163*ql[0]*D[11]-0.3977475644174328*qc[0]*D[11]-0.3314563036811939*qr[6]*D[10]-0.3314563036811939*ql[6]*D[10]-0.6629126073623879*qc[6]*D[10]+0.3444594950788841*qr[3]*D[10]-0.3444594950788841*ql[3]*D[10]-0.1913663861549356*qr[7]*D[9]+0.1913663861549356*ql[7]*D[9]+0.1988737822087163*qr[5]*D[9]+0.1988737822087163*ql[5]*D[9]-0.3977475644174328*qc[5]*D[9]-0.1913663861549356*qr[6]*D[8]+0.1913663861549356*ql[6]*D[8]+0.1988737822087163*qr[3]*D[8]+0.1988737822087163*ql[3]*D[8]-0.3977475644174328*qc[3]*D[8]);
  out[4] += J*((-0.2679129406169099*qr[6]*D[15])+0.2679129406169099*ql[6]*D[15]+0.3314563036811939*qr[3]*D[15]+0.3314563036811939*ql[3]*D[15]-0.6629126073623879*qc[3]*D[15]-0.2679129406169099*qr[7]*D[14]+0.2679129406169099*ql[7]*D[14]+0.3314563036811939*qr[5]*D[14]+0.3314563036811939*ql[5]*D[14]-0.6629126073623879*qc[5]*D[14]-0.1546796083845572*qr[6]*D[13]-0.1546796083845572*ql[6]*D[13]-1.016465997955661*qc[6]*D[13]+0.1913663861549356*qr[3]*D[13]-0.1913663861549356*ql[3]*D[13]-0.2679129406169099*qr[2]*D[12]+0.2679129406169099*ql[2]*D[12]+0.3314563036811939*qr[0]*D[12]+0.3314563036811939*ql[0]*D[12]-0.6629126073623879*qc[0]*D[12]-0.1546796083845572*qr[7]*D[11]-0.1546796083845572*ql[7]*D[11]-1.016465997955661*qc[7]*D[11]+0.1913663861549356*qr[5]*D[11]-0.1913663861549356*ql[5]*D[11]-0.2679129406169099*qr[4]*D[10]+0.2679129406169099*ql[4]*D[10]+0.3314563036811939*qr[1]*D[10]+0.3314563036811939*ql[1]*D[10]-0.6629126073623879*qc[1]*D[10]-0.1546796083845572*qr[2]*D[9]-0.1546796083845572*ql[2]*D[9]-1.016465997955661*qc[2]*D[9]+0.1913663861549356*qr[0]*D[9]-0.1913663861549356*ql[0]*D[9]-0.1546796083845572*qr[4]*D[8]-0.1546796083845572*ql[4]*D[8]-1.016465997955661*qc[4]*D[8]+0.1913663861549356*qr[1]*D[8]-0.1913663861549356*ql[1]*D[8]);
  out[5] += J*((-0.3314563036811939*qr[2]*D[15])-0.3314563036811939*ql[2]*D[15]-0.6629126073623879*qc[2]*D[15]+0.3444594950788841*qr[0]*D[15]-0.3444594950788841*ql[0]*D[15]-0.3314563036811939*qr[4]*D[14]-0.3314563036811939*ql[4]*D[14]-0.6629126073623879*qc[4]*D[14]+0.3444594950788841*qr[1]*D[14]-0.3444594950788841*ql[1]*D[14]-0.1913663861549356*qr[2]*D[13]+0.1913663861549356*ql[2]*D[13]+0.1988737822087163*qr[0]*D[13]+0.1988737822087163*ql[0]*D[13]-0.3977475644174328*qc[0]*D[13]-0.3314563036811939*qr[6]*D[12]-0.3314563036811939*ql[6]*D[12]-0.6629126073623879*qc[6]*D[12]+0.3444594950788841*qr[3]*D[12]-0.3444594950788841*ql[3]*D[12]-0.1913663861549356*qr[4]*D[11]+0.1913663861549356*ql[4]*D[11]+0.1988737822087163*qr[1]*D[11]+0.1988737822087163*ql[1]*D[11]-0.3977475644174328*qc[1]*D[11]-0.3314563036811939*qr[7]*D[10]-0.3314563036811939*ql[7]*D[10]-0.6629126073623879*qc[7]*D[10]+0.3444594950788841*qr[5]*D[10]-0.3444594950788841*ql[5]*D[10]-0.1913663861549356*qr[6]*D[9]+0.1913663861549356*ql[6]*D[9]+0.1988737822087163*qr[3]*D[9]+0.1988737822087163*ql[3]*D[9]-0.3977475644174328*qc[3]*D[9]-0.1913663861549356*qr[7]*D[8]+0.1913663861549356*ql[7]*D[8]+0.1988737822087163*qr[5]*D[8]+0.1988737822087163*ql[5]*D[8]-0.3977475644174328*qc[5]*D[8]);
  out[6] += J*((-0.2679129406169099*qr[4]*D[15])+0.2679129406169099*ql[4]*D[15]+0.3314563036811939*qr[1]*D[15]+0.3314563036811939*ql[1]*D[15]-0.6629126073623879*qc[1]*D[15]-0.2679129406169099*qr[2]*D[14]+0.2679129406169099*ql[2]*D[14]+0.3314563036811939*qr[0]*D[14]+0.3314563036811939*ql[0]*D[14]-0.6629126073623879*qc[0]*D[14]-0.1546796083845572*qr[4]*D[13]-0.1546796083845572*ql[4]*D[13]-1.016465997955661*qc[4]*D[13]+0.1913663861549356*qr[1]*D[13]-0.1913663861549356*ql[1]*D[13]-0.2679129406169099*qr[7]*D[12]+0.2679129406169099*ql[7]*D[12]+0.3314563036811939*qr[5]*D[12]+0.3314563036811939*ql[5]*D[12]-0.6629126073623879*qc[5]*D[12]-0.1546796083845572*qr[2]*D[11]-0.1546796083845572*ql[2]*D[11]-1.016465997955661*qc[2]*D[11]+0.1913663861549356*qr[0]*D[11]-0.1913663861549356*ql[0]*D[11]-0.2679129406169099*qr[6]*D[10]+0.2679129406169099*ql[6]*D[10]+0.3314563036811939*qr[3]*D[10]+0.3314563036811939*ql[3]*D[10]-0.6629126073623879*qc[3]*D[10]-0.1546796083845572*qr[7]*D[9]-0.1546796083845572*ql[7]*D[9]-1.016465997955661*qc[7]*D[9]+0.1913663861549356*qr[5]*D[9]-0.1913663861549356*ql[5]*D[9]-0.1546796083845572*qr[6]*D[8]-0.1546796083845572*ql[6]*D[8]-1.016465997955661*qc[6]*D[8]+0.1913663861549356*qr[3]*D[8]-0.1913663861549356*ql[3]*D[8]);
  out[7] += J*((-0.2679129406169099*qr[2]*D[15])+0.2679129406169099*ql[2]*D[15]+0.3314563036811939*qr[0]*D[15]+0.3314563036811939*ql[0]*D[15]-0.6629126073623879*qc[0]*D[15]-0.2679129406169099*qr[4]*D[14]+0.2679129406169099*ql[4]*D[14]+0.3314563036811939*qr[1]*D[14]+0.3314563036811939*ql[1]*D[14]-0.6629126073623879*qc[1]*D[14]-0.1546796083845572*qr[2]*D[13]-0.1546796083845572*ql[2]*D[13]-1.016465997955661*qc[2]*D[13]+0.1913663861549356*qr[0]*D[13]-0.1913663861549356*ql[0]*D[13]-0.2679129406169099*qr[6]*D[12]+0.2679129406169099*ql[6]*D[12]+0.3314563036811939*qr[3]*D[12]+0.3314563036811939*ql[3]*D[12]-0.6629126073623879*qc[3]*D[12]-0.1546796083845572*qr[4]*D[11]-0.1546796083845572*ql[4]*D[11]-1.016465997955661*qc[4]*D[11]+0.1913663861549356*qr[1]*D[11]-0.1913663861549356*ql[1]*D[11]-0.2679129406169099*qr[7]*D[10]+0.2679129406169099*ql[7]*D[10]+0.3314563036811939*qr[5]*D[10]+0.3314563036811939*ql[5]*D[10]-0.6629126073623879*qc[5]*D[10]-0.1546796083845572*qr[6]*D[9]-0.1546796083845572*ql[6]*D[9]-1.016465997955661*qc[6]*D[9]+0.1913663861549356*qr[3]*D[9]-0.1913663861549356*ql[3]*D[9]-0.1546796083845572*qr[7]*D[8]-0.1546796083845572*ql[7]*D[8]-1.016465997955661*qc[7]*D[8]+0.1913663861549356*qr[5]*D[8]-0.1913663861549356*ql[5]*D[8]);
}
