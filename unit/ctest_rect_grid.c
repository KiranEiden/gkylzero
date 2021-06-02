#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_rect_grid.h>

void test_grid_2d()
{
  double lower[] = {1.0, 1.0}, upper[] = {2.5, 5.0};
  int cells[] = {20, 20};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  TEST_CHECK( grid.ndim == 2 );
  for (int i=0; i<grid.ndim; ++i) {
    TEST_CHECK( grid.lower[i] == lower[i] );
    TEST_CHECK( grid.upper[i] == upper[i] );
    TEST_CHECK( grid.cells[i] == cells[i] );
    TEST_CHECK( grid.dx[i] == (upper[i]-lower[i])/cells[i] );
  }
  TEST_CHECK( grid.cellVolume == 0.075*0.2 );

  int idx[2];
  double xc[2];
  for (int i=0; i<grid.cells[0]; ++i)
    for (int j=0; j<grid.cells[1]; ++j) {
      idx[0] = i; idx[1] = j;
      gkyl_rect_grid_cell_center(&grid, idx, xc);

      TEST_CHECK( xc[0] == 1.0 + (i+0.5)*grid.dx[0] );
      TEST_CHECK( xc[1] == 1.0 + (j+0.5)*grid.dx[1] );
    }
}

// CUDA specific tests
#ifdef GKYL_HAVE_CUDA

int cu_rect_grid_test(const struct gkyl_rect_grid *rng);

void test_cu_grid_2d()
{
  double lower[] = {1.0, 1.0}, upper[] = {2.5, 5.0};
  int cells[] = {20, 20};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  // clone on device
  struct gkyl_rect_grid *grid_cu = gkyl_rect_grid_clone_on_cu_dev(&grid);

  int nfail = cu_rect_grid_test(grid_cu);
  TEST_CHECK( nfail == 0 );
  
  gkyl_cu_free(grid_cu);
}
#endif

TEST_LIST = {
  { "grid_2d", test_grid_2d },
#ifdef GKYL_HAVE_CUDA
  { "cu_grid_2d", test_cu_grid_2d },
#endif  
  { NULL, NULL },
};
