#pragma once

#include <gkyl_array.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>
#include <gkyl_ref_count.h>


#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

#ifdef GKYL_HAVE_CUDA
//#define CUBLAS_H_
//#include <cuda_runtime.h>
#if !defined(CUBLAS_V2_H_)
#include <cublas_v2.h>
#endif

/**
 * Computes: alpha*matrix_multiplication(A,B) + Beta*C = C 
 * This is done using the cublas_v2 library. The function here, cu_mat_mm_array is designed specifically so
 * A is type gkyl_mat_trans, and B/C are gkyl_arrays. This is purposefully done for calculations with
 * small A but and very large B, C where we don't wish to translate B, C to mat/nmat types. 
 * 
 * Such calculations use this function like the conversion from nodal to modal representation of phase space
 * quantities.
 * @param cuh cublasHandle_t object
 * @param alpha Coefficient infront of A*B
 * @param beta Coefficient infron of C
 * @param transa Whether or not to transpose A
 * @param A gkyl_mat matrix for computing A*B = C
 * @param transb Whether or not to transpose B
 * @param B gkyl_array matrix for computing A*B = C
 * @param C gkyl_array matrix for computing A*B = C
*/
void cu_mat_mm_array(cublasHandle_t cuh, double alpha, double beta, enum gkyl_mat_trans transa, struct gkyl_mat *A, 
  enum gkyl_mat_trans transb, struct gkyl_array *B, struct gkyl_array *C);
  
#endif

struct gkyl_phase_nodal_to_modal_mem {

  // info for alpha*matrix_multiplication(A,B) + Beta*C = C 
  // using cu_mat_mm_array
  double alpha;
  double beta;
  enum gkyl_mat_trans transa;
  enum gkyl_mat_trans transb;
  //struct gkyl_mat *A_ho;
  //struct gkyl_mat *A_cu;

#ifdef GKYL_HAVE_CUDA
  cublasHandle_t cuh; // cublas handle
#endif  
};