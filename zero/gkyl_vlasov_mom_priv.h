#pragma once

// Private header, not for direct use in user code

#include <gkyl_vlasov_mom_kernels.h>

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

typedef void (*vlasov_momf_t)(const double *xc, const double *dx,
  const int *idx, const double *fIn, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vlasov_momf_t kernels[3]; } gkyl_mom_kern_list;

//
// Serendipity basis kernels
//

// M0 kernel list
GKYL_CU_D
static const gkyl_mom_kern_list ser_m0_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M0_1x1v_ser_p1, vlasov_M0_1x1v_ser_p2 }, // 0
  { NULL, vlasov_M0_1x2v_ser_p1, vlasov_M0_1x2v_ser_p2 }, // 1
  { NULL, vlasov_M0_1x3v_ser_p1, vlasov_M0_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M0_2x2v_ser_p1, vlasov_M0_2x2v_ser_p2 }, // 3
  { NULL, vlasov_M0_2x3v_ser_p1, NULL                  }, // 4
  // 3x kernels
  { NULL, vlasov_M0_3x3v_ser_p1, NULL                  }, // 5
};

// M1i kernel list
GKYL_CU_D
static const gkyl_mom_kern_list ser_m1i_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M1i_1x1v_ser_p1, vlasov_M1i_1x1v_ser_p2 }, // 0
  { NULL, vlasov_M1i_1x2v_ser_p1, vlasov_M1i_1x2v_ser_p2 }, // 1
  { NULL, vlasov_M1i_1x3v_ser_p1, vlasov_M1i_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M1i_2x2v_ser_p1, vlasov_M1i_2x2v_ser_p2 }, // 3
  { NULL, vlasov_M1i_2x3v_ser_p1, NULL                   }, // 4
  // 3x kernels
  { NULL, vlasov_M1i_3x3v_ser_p1, NULL                   }, // 5
};

// M2 kernel list
GKYL_CU_D
static const gkyl_mom_kern_list ser_m2_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M2_1x1v_ser_p1, vlasov_M2_1x1v_ser_p2 }, // 0
  { NULL, vlasov_M2_1x2v_ser_p1, vlasov_M2_1x2v_ser_p2 }, // 1
  { NULL, vlasov_M2_1x3v_ser_p1, vlasov_M2_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M2_2x2v_ser_p1, vlasov_M2_2x2v_ser_p2 }, // 3
  { NULL, vlasov_M2_2x3v_ser_p1, NULL                  }, // 4
  // 3x kernels
  { NULL, vlasov_M2_3x3v_ser_p1, NULL                  }, // 5
};

// M2ij kernel list
GKYL_CU_D
static const gkyl_mom_kern_list ser_m2ij_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M2ij_1x1v_ser_p1, vlasov_M2ij_1x1v_ser_p2 }, // 0
  { NULL, vlasov_M2ij_1x2v_ser_p1, vlasov_M2ij_1x2v_ser_p2 }, // 1
  { NULL, vlasov_M2ij_1x3v_ser_p1, vlasov_M2ij_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M2ij_2x2v_ser_p1, vlasov_M2ij_2x2v_ser_p2 }, // 3
  { NULL, vlasov_M2ij_2x3v_ser_p1, NULL                    }, // 4
  // 3x kernels
  { NULL, vlasov_M2ij_3x3v_ser_p1, NULL                    }, // 5
};

// M3i kernel list
GKYL_CU_D
static const gkyl_mom_kern_list ser_m3i_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M3i_1x1v_ser_p1, vlasov_M3i_1x1v_ser_p2 }, // 0
  { NULL, vlasov_M3i_1x2v_ser_p1, vlasov_M3i_1x2v_ser_p2 }, // 1
  { NULL, vlasov_M3i_1x3v_ser_p1, vlasov_M3i_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M3i_2x2v_ser_p1, vlasov_M3i_2x2v_ser_p2 }, // 3
  { NULL, vlasov_M3i_2x3v_ser_p1, NULL                   }, // 4
  // 3x kernels
  { NULL, vlasov_M3i_3x3v_ser_p1, NULL                   }, // 5
};

// M3ijk kernel list
GKYL_CU_D
static const gkyl_mom_kern_list ser_m3ijk_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M3ijk_1x1v_ser_p1, vlasov_M3ijk_1x1v_ser_p2 }, // 0
  { NULL, vlasov_M3ijk_1x2v_ser_p1, vlasov_M3ijk_1x2v_ser_p2 }, // 1
  { NULL, vlasov_M3ijk_1x3v_ser_p1, vlasov_M3ijk_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M3ijk_2x2v_ser_p1, vlasov_M3ijk_2x2v_ser_p2 }, // 3
  { NULL, vlasov_M3ijk_2x3v_ser_p1, NULL                     }, // 4
  // 3x kernels
  { NULL, vlasov_M3ijk_3x3v_ser_p1, NULL                     }, // 5
};

//
// Tensor-product basis kernels
//

// M0 kernel list
GKYL_CU_D
static const gkyl_mom_kern_list ten_m0_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M0_1x1v_ser_p1, vlasov_M0_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_M0_1x2v_ser_p1, vlasov_M0_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_M0_1x3v_ser_p1, vlasov_M0_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M0_2x2v_ser_p1, vlasov_M0_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_M0_2x3v_ser_p1, NULL                  }, // 4
  // 3x kernels
  { NULL, vlasov_M0_3x3v_ser_p1, NULL                  }, // 5
};

// M1i kernel list
GKYL_CU_D
static const gkyl_mom_kern_list ten_m1i_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M1i_1x1v_ser_p1, vlasov_M1i_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_M1i_1x2v_ser_p1, vlasov_M1i_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_M1i_1x3v_ser_p1, vlasov_M1i_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M1i_2x2v_ser_p1, vlasov_M1i_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_M1i_2x3v_ser_p1, NULL                   }, // 4
  // 3x kernels
  { NULL, vlasov_M1i_3x3v_ser_p1, NULL                   }, // 5
};

// M2 kernel list
GKYL_CU_D
static const gkyl_mom_kern_list ten_m2_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M2_1x1v_ser_p1, vlasov_M2_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_M2_1x2v_ser_p1, vlasov_M2_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_M2_1x3v_ser_p1, vlasov_M2_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M2_2x2v_ser_p1, vlasov_M2_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_M2_2x3v_ser_p1, NULL                  }, // 4
  // 3x kernels
  { NULL, vlasov_M2_3x3v_ser_p1, NULL                  }, // 5
};

// M2ij kernel list
GKYL_CU_D
static const gkyl_mom_kern_list ten_m2ij_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M2ij_1x1v_ser_p1, vlasov_M2ij_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_M2ij_1x2v_ser_p1, vlasov_M2ij_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_M2ij_1x3v_ser_p1, vlasov_M2ij_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M2ij_2x2v_ser_p1, vlasov_M2ij_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_M2ij_2x3v_ser_p1, NULL                    }, // 4
  // 3x kernels
  { NULL, vlasov_M2ij_3x3v_ser_p1, NULL                    }, // 5
};

// M3i kernel list
GKYL_CU_D
static const gkyl_mom_kern_list ten_m3i_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M3i_1x1v_ser_p1, vlasov_M3i_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_M3i_1x2v_ser_p1, vlasov_M3i_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_M3i_1x3v_ser_p1, vlasov_M3i_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M3i_2x2v_ser_p1, vlasov_M3i_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_M3i_2x3v_ser_p1, NULL                   }, // 4
  // 3x kernels
  { NULL, vlasov_M3i_3x3v_ser_p1, NULL                   }, // 5
};

// M3ijk kernel list
GKYL_CU_D
static const gkyl_mom_kern_list ten_m3ijk_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M3ijk_1x1v_ser_p1, vlasov_M3ijk_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_M3ijk_1x2v_ser_p1, vlasov_M3ijk_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_M3ijk_1x3v_ser_p1, vlasov_M3ijk_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M3ijk_2x2v_ser_p1, vlasov_M3ijk_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_M3ijk_2x3v_ser_p1, NULL                     }, // 4
  // 3x kernels
  { NULL, vlasov_M3ijk_3x3v_ser_p1, NULL                     }, // 5
};

struct vlasov_mom_type {
  struct gkyl_mom_type momt;
  int cdim; // config-space dim
  int pdim; // phase-space dim
  int poly_order; // polynomal order
  int num_config; // number of basis functions in config-space
  int num_phase; // number of basis functions in phase-space
  int num_mom; // number of components in moment
  vlasov_momf_t kernel; // moment calculation kernel
  struct gkyl_ref_count ref_count; // reference count
};

GKYL_CU_D
static void
kernel(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out)
{
  struct vlasov_mom_type *vlasov_mom = container_of(momt, struct vlasov_mom_type, momt);

  return vlasov_mom->kernel(xc, dx, idx, f, out);
}
