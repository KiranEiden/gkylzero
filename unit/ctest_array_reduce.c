#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <time.h>

void test_dummy()
{
}

// CUDA specific tests
#ifdef GKYL_HAVE_CUDA

void test_cu_array_reduce_max()
{

  unsigned long numComp = 3, numCells = 10;
  // Create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  // Initialize data
  double *a1_d = a1->data;
  for (unsigned i=0; i<numCells; ++i) {
    for (unsigned k=0; k<numComp; ++k) {
      a1_d[i*numComp+k] = (double)i+(double)k*0.10;
  }}
  gkyl_array_copy(a1_cu, a1);
  // create a device and host arrays to store reduction
  double* a1max = gkyl_malloc(numComp*sizeof(double));
  double* a1max_cu = (double*) gkyl_cu_malloc(numComp*sizeof(double));

  // Component-wise reduce array.
  gkyl_array_reduce(a1max_cu, a1_cu, GKYL_MAX);

  // Copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k=0; k<numComp; ++k) {
    TEST_CHECK( gkyl_compare(a1max[k], (double)(numCells-1)+(double)k*0.10, 1e-14) );
  }

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);

}

void test_cu_array_reduce_range_1d_max()
{
  unsigned long numComp = 1, numCells = 10;
  // create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  // initialize data
  double *a1_d = a1->data;
  for (unsigned i=0; i<numCells; ++i) {
    for (unsigned k=0; k<numComp; ++k) {
      a1_d[i*numComp+k] = (double)i+(double)k*0.10;
  }}
  gkyl_array_copy(a1_cu, a1);
  // create a device and host arrays to store reduction
  double* a1max = gkyl_malloc(numComp*sizeof(double));
  double* a1max_cu = (double*) gkyl_cu_malloc(numComp*sizeof(double));

  int lower[] = {1}, upper[] = {numCells};
  int sublower[] = {3}, subupper[] = {7};
  struct gkyl_range range, subrange;
  gkyl_range_init(&range, 1, lower, upper);
  gkyl_sub_range_init(&subrange, &range, sublower, subupper);

  // Component-wise reduce array on range
  gkyl_array_reduce_range(a1max_cu, a1_cu, GKYL_MAX, range);
  // Copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k=0; k<numComp; ++k) {
    TEST_CHECK( gkyl_compare(a1max[k], (double)(upper[0]-1)+(double)k*0.10, 1e-14) );
  }

  // Component-wise reduce array on subrange
  gkyl_cu_memset(a1max_cu, 0, sizeof(double)*numComp);
  gkyl_array_reduce_range(a1max_cu, a1_cu, GKYL_MAX, subrange);

  // Copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k=0; k<numComp; ++k) {
    TEST_CHECK( gkyl_compare(a1max[k], (double)(subupper[0]-1)+(double)k*0.10, 1e-14) );
  }

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);
}

void test_cu_array_reduce_range_2d_max()
{
  int cells[] = {8, 10};
  int ghost[] = {1, 0};
  double lower[] = {0., -1.};
  double upper[] = {1., 1.};

  struct gkyl_rect_grid grid;
  struct gkyl_range range, range_ext;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);
  gkyl_create_grid_ranges(&grid, ghost, &range_ext, &range);

  unsigned long numComp = 1;
  unsigned long numCells = range_ext.volume;
  // create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *b1 = gkyl_array_new(GKYL_DOUBLE, 1, numCells);
  struct gkyl_array *b1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, numCells);
  // initialize data
  double *a1_d = a1->data;
  double *b1_d = b1->data;
  for (unsigned i=0; i<numCells; ++i) {
    b1_d[i] = (double) i+1000;
    for (unsigned k=0; k<numComp; ++k) {
      a1_d[i*numComp+k] = (double)numCells-1-i+(double)k*0.10;
  }}
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(b1_cu, b1);
  // create a device and host arrays to store reduction
  double* a1max = gkyl_malloc(numComp*sizeof(double));
  double* a1max_correct = gkyl_malloc(numComp*sizeof(double));
  double* a1max_cu = (double*) gkyl_cu_malloc(numComp*sizeof(double));
  double* b1max = gkyl_malloc(1*sizeof(double));
  double* b1max_cu = (double*) gkyl_cu_malloc(1*sizeof(double));
  gkyl_cu_memset(a1max_cu, 0, sizeof(double)*numComp);
  gkyl_cu_memset(b1max_cu, 0, sizeof(double)*1);

  // Component-wise reduce array on range
  gkyl_array_reduce_range(b1max_cu, b1_cu, GKYL_MAX, range_ext);
  // Copy to host and check values.
  gkyl_cu_memcpy(b1max, b1max_cu, sizeof(double), GKYL_CU_MEMCPY_D2H);
  TEST_CHECK( gkyl_compare(b1max[0], (double)(1000+numCells-1), 1e-14) );

  // Component-wise reduce array on subrange
  gkyl_cu_memset(a1max_cu, 0, sizeof(double)*numComp);
  gkyl_array_reduce_range(a1max_cu, a1_cu, GKYL_MAX, range);

  gkyl_array_reduce_range(a1max_correct, a1, GKYL_MAX, range);

  // Copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k=0; k<numComp; ++k) {
    TEST_CHECK( gkyl_compare(a1max[k], a1max_correct[k], 1e-14) );
  }

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);
}

void
test_cu_array_reduce_range_max_timer(int NX, int NY, int VX, int VY)
{
  int cells[] = {NX, NY, VX, VY};
  int ghost[] = {1, 1, 0, 0};
  double lower[] = {0.0, 0.0, 0.0, 0.0};
  double upper[] = {1.0, 1.0, 1.0, 1.0};

  int ndim = 4;
  struct gkyl_rect_grid grid;
  struct gkyl_range range, range_ext;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  gkyl_create_grid_ranges(&grid, ghost, &range_ext, &range);

  unsigned long numComp = 1;
  unsigned long numCells = range_ext.volume;
  // create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *b1 = gkyl_array_new(GKYL_DOUBLE, 1, numCells);
  struct gkyl_array *b1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, numCells);
  // initialize data
  double *a1_d = a1->data;
  double *b1_d = b1->data;
  for (unsigned i=0; i<numCells; ++i) {
    b1_d[i] = (double) i+1000;
    for (unsigned k=0; k<numComp; ++k) {
      a1_d[i*numComp+k] = (double)numCells-1-i+(double)k*0.10;
  }}
  gkyl_array_copy(a1_cu, a1);
  gkyl_array_copy(b1_cu, b1);
  // create a device and host arrays to store reduction
  double* a1max = gkyl_malloc(numComp*sizeof(double));
  double* a1max_correct = gkyl_malloc(numComp*sizeof(double));
  double* a1max_cu = (double*) gkyl_cu_malloc(numComp*sizeof(double));
  double* b1max = gkyl_malloc(1*sizeof(double));
  double* b1max_cu = (double*) gkyl_cu_malloc(1*sizeof(double));
  gkyl_cu_memset(a1max_cu, 0, sizeof(double)*numComp);
  gkyl_cu_memset(b1max_cu, 0, sizeof(double)*1);

  // Component-wise reduce array on range
  gkyl_array_reduce_range(b1max_cu, b1_cu, GKYL_MAX, range_ext);
  // Copy to host and check values.
  gkyl_cu_memcpy(b1max, b1max_cu, sizeof(double), GKYL_CU_MEMCPY_D2H);
  TEST_CHECK( gkyl_compare(b1max[0], (double)(1000+numCells-1), 1e-14) );

  // Component-wise reduce array on subrange
  gkyl_cu_memset(a1max_cu, 0, sizeof(double)*numComp);

  struct timespec tm = gkyl_wall_clock();

  for (int i=0; i<100; ++i)
    gkyl_array_reduce_range(a1max_cu, a1_cu, GKYL_MAX, range);
  double red_tm = gkyl_time_diff_now_sec(tm);

  printf("100 reductions on (%d,%d,%d,%d) took %g sec\n", NX, NY, VX, VY,
    red_tm);

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);
  gkyl_array_release(b1);
  gkyl_array_release(b1_cu);
}

void
test_cu_array_reduce_range_max_timer_32x32x40x40()
{
  test_cu_array_reduce_range_max_timer(32, 32, 40, 40);
}

void
test_cu_array_reduce_range_max_timer_32x32x32x32()
{
  test_cu_array_reduce_range_max_timer(32, 32, 32, 32);
}

#endif

TEST_LIST = {
  { "dummy", test_dummy },
#ifdef GKYL_HAVE_CUDA
  { "cu_array_reduce_max", test_cu_array_reduce_max },
  { "cu_array_reduce_range_1d_max", test_cu_array_reduce_range_1d_max },
  { "cu_array_reduce_range_2d_max", test_cu_array_reduce_range_2d_max },
  { "cu_array_reduce_range_max_timer_32x32x40x40", test_cu_array_reduce_range_max_timer_32x32x40x40  },
  { "cu_array_reduce_range_max_timer_32x32x32x32", test_cu_array_reduce_range_max_timer_32x32x32x32  },  
#endif
  { NULL, NULL },
};
