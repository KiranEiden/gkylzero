#include <acutest.h>

#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

void evalFunc(double t, const double *xn, double* restrict fout)
{
  double x = xn[0];
  fout[0] = x*x;
}

void
test_1()
{
  int polyOrder = 1;
  double lower[] = {-2.0}, upper[] = {2.0};
  int cells[] = {2};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, polyOrder);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    polyOrder+1, 1, evalFunc);

  // create array range: no ghost-cells in velocity space
  struct gkyl_range arr_range;
  gkyl_range_init_from_shape(&arr_range, 1, cells);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(sizeof(double)*basis.numBasis,
    arr_range.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);

  for (unsigned i=0; i<distf->size; ++i) {
    const double *df = gkyl_array_fetch(distf, i);
    printf("%g %g\n", df[0], df[1]);
  }
  printf("\n");

}

TEST_LIST = {
  { "test_1", test_1 },
  { NULL, NULL },
};
