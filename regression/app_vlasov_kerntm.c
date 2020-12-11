#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

#define SHOW_TIME(msg, tdiff) printf("%s %g s\n", msg, 1.0*(tdiff)/CLOCKS_PER_SEC)

void
init_phase_ranges(int cdim, int pdim, const int *cells,
  struct gkyl_range *phase_local, struct gkyl_range *phase_local_ext)
{
  int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  
  for (unsigned i=0; i<cdim; ++i) {
    lower_ext[i] = -1;
    upper_ext[i] = cells[i];

    lower[i] = 0;
    upper[i] = cells[i]-1;
  }
  for (unsigned i=cdim; i<pdim; ++i) {
    lower_ext[i] = 0;
    upper_ext[i] = cells[i]-1;

    lower[i] = 0;
    upper[i] = cells[i]-1;
  }
  gkyl_range_init(phase_local_ext, pdim, lower_ext, upper_ext);
  gkyl_sub_range_init(phase_local, phase_local_ext, lower, upper);
}

void
init_conf_ranges(int cdim, const int *cells,
  struct gkyl_range *conf_local, struct gkyl_range *conf_local_ext)
{
  int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  
  for (unsigned i=0; i<cdim; ++i) {
    lower_ext[i] = -1;
    upper_ext[i] = cells[i];

    lower[i] = 0;
    upper[i] = cells[i]-1;
  }
  gkyl_range_init(conf_local_ext, cdim, lower_ext, upper_ext);
  gkyl_sub_range_init(conf_local, conf_local_ext, lower, upper);
}

struct gkyl_array*
mkarr(long numComp, long size)
{
  return gkyl_array_new(sizeof(double)*numComp, size);
}

struct vlasov_app {
    int cdim, vdim, polyOrder;
    int ccells[3], vcells[3];
    int nloop;
};

struct vlasov_app
parse_args(int argc, char **argv)
{
  // replace with something read from CLI or a config file
  struct vlasov_app app = {
    .cdim = 2, .vdim = 2, .polyOrder = 2,
    .ccells = {8, 8, 8}, .vcells = {16, 16, 16},
    .nloop = 100
  };

  // struct vlasov_app app = {
  //   .cdim = 2, .vdim = 2, .polyOrder = 2,
  //   .ccells = {2, 2}, .vcells = {4, 4},
  //   .nloop = 1
  // };

  return app;
}

static inline void
copy_int_arr(int n, const int *restrict inp, int *restrict out)
{
  for (unsigned i=0; i<n; ++i) out[i] = inp[i];
}

void
show_params(struct vlasov_app app)
{
  int cdim = app.cdim, vdim = app.vdim;
  int polyOrder = app.polyOrder;  
  
  printf("Kernel timer run with:\n");
  printf("cdim = %d; vdim = %d; polyOrder = %d\n", cdim, vdim, polyOrder);
  printf("nloop = %d\n", app.nloop);
  
  if (cdim == 1)
    printf("NX = %d\n", app.ccells[0]);
  if (cdim == 2)
    printf("NX = %d; NY = %d\n", app.ccells[0], app.ccells[1]);
  if (cdim == 3)
    printf("NX = %d; NY = %d; NZ = %d\n", app.ccells[0], app.ccells[1], app.ccells[2]);

  if (vdim == 1)
    printf("VX = %d\n", app.vcells[0]);
  if (vdim == 2)
    printf("VX = %d; VY = %d\n", app.vcells[0], app.vcells[1]);
  if (vdim == 3)
    printf("VX = %d; VY = %d; VZ = %d\n", app.vcells[0], app.vcells[1], app.vcells[2]);
}

int
main(int argc, char **argv)
{
  struct vlasov_app app = parse_args(argc, argv);
  show_params(app);

  int cdim = app.cdim, vdim = app.vdim, pdim = cdim+vdim;
  int polyOrder = app.polyOrder;
  int nloop = app.nloop;

  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM];
  for (unsigned d=0; d<cdim; ++d) {
    lower[d] = -1.0; upper[d] = 1.0;
    cells[d] = app.ccells[d];
  }
  for (unsigned d=cdim; d<pdim; ++d) {
    lower[d] = -6.0; upper[d] = 6.0;
    cells[d] = app.vcells[d-cdim];
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);  

  // basis functions
  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, pdim, polyOrder);
  gkyl_cart_modal_serendip(&confBasis, cdim, polyOrder);

  struct gkyl_range phase_local, phase_local_ext, conf_local, conf_local_ext;
  init_phase_ranges(cdim, pdim, cells, &phase_local, &phase_local_ext);
  init_conf_ranges(cdim, cells, &conf_local, &conf_local_ext);

  // create distribution function, fields
  struct gkyl_array *fIn = mkarr(basis.numBasis, phase_local_ext.volume);
  struct gkyl_array *fOut = mkarr(basis.numBasis, phase_local_ext.volume);
  struct gkyl_array *em = mkarr(confBasis.numBasis*8, conf_local_ext.volume);

  gkyl_array_clear(fIn, 1.0); gkyl_array_clear(fOut, 1.0);
  gkyl_array_clear(em, 2.0);

  // vlasov equation object
  struct gkyl_dg_eqn *vlasov_eqn = gkyl_dg_vlasov_new(&confBasis, &basis, &conf_local);

  clock_t tstart, tend;
  double xc[GKYL_MAX_DIM];

  long nvol = 0;

  // volume kernels
  printf("\n*** Timing volume kernel ....\n\n");  
  tstart = clock();
  for (unsigned n=0; n<nloop; ++n) {
    
    gkyl_vlasov_set_qmem(vlasov_eqn, em);

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &phase_local);    
    while(gkyl_range_iter_next(&iter)) {
      gkyl_rect_grid_cell_center(&grid, iter.idx, xc);
      
      long fidx = gkyl_range_idx(&phase_local, iter.idx);
      
      vlasov_eqn->vol_term(vlasov_eqn, xc, grid.dx, iter.idx,
        gkyl_array_fetch(fIn, fidx), gkyl_array_fetch(fOut, fidx)
      );
      nvol += 1;
    }
  }
  tend = clock();
  SHOW_TIME("Volume kernel took ", tend-tstart);
  printf("(Total calls = %ld. Time per-call %g)\n", nvol,
    1.0*(tend-tstart)/CLOCKS_PER_SEC/nvol);

  double totVolTm = 1.0*(tend-tstart)/CLOCKS_PER_SEC;

  // surface kernels
  printf("\n*** Timing surface kernels ....\n\n");

  int zero_flux_dir[6];
  for (unsigned i=0; i<cdim; ++i)
    zero_flux_dir[i] = 0; // ghost cells in config directions
  for (unsigned i=cdim; i<pdim; ++i)
    zero_flux_dir[i] = 1; // no ghost cells in velocity directions

  long nsurf[] = { 0, 0, 0, 0, 0, 0 };
  clock_t tmsurf[6];

  for (unsigned n=0; n<nloop; ++n) {
    gkyl_vlasov_set_qmem(vlasov_eqn, em);

    for (unsigned dir=0; dir<pdim; ++dir) {
      tstart = clock();

      int dirLoIdx = phase_local.lower[dir];
      int dirUpIdx = phase_local.upper[dir]+1; // one more edge than cells

      if (zero_flux_dir[dir]) {
        dirLoIdx = dirLoIdx+1;
        dirUpIdx = dirUpIdx-1;
      }

      struct gkyl_range perp_range;
      gkyl_range_shorten(&perp_range, &phase_local, dir, 1);

      struct gkyl_range_iter iter;
      gkyl_range_iter_init(&iter, &perp_range);

      int idxm[GKYL_MAX_DIM], idxp[GKYL_MAX_DIM];
      double xcm[GKYL_MAX_DIM], xcp[GKYL_MAX_DIM];

      while(gkyl_range_iter_next(&iter)) {
        copy_int_arr(pdim, iter.idx, idxm);
        copy_int_arr(pdim, iter.idx, idxp);

        for (int i=dirLoIdx; i<=dirUpIdx; ++i) { // note dirUpIdx is inclusive
          idxm[dir] = i-1; idxp[dir] = i;
          
          gkyl_rect_grid_cell_center(&grid, idxm, xcm);
          gkyl_rect_grid_cell_center(&grid, idxp, xcp);

          long linIdxm = gkyl_range_idx(&phase_local, idxm);
          long linIdxp = gkyl_range_idx(&phase_local, idxp);

          vlasov_eqn->surf_term(vlasov_eqn, dir,
            xcm, xcp, grid.dx, grid.dx, 1.0, idxm, idxp,
            gkyl_array_fetch(fIn, linIdxm), gkyl_array_fetch(fIn, linIdxp),
            gkyl_array_fetch(fOut, linIdxm), gkyl_array_fetch(fOut, linIdxp)
          );
          nsurf[dir] += 1;
        }
      }
      
      tend = clock();
      tmsurf[dir] = tend-tstart;
    }
  }

  double totSurfTm = 0.0;
  for (unsigned dir=0; dir<pdim; ++dir) {
    printf("Surface update in dir %d took %g sec\n",
      dir, 1.0*tmsurf[dir]/CLOCKS_PER_SEC);
    printf(" (Total calls = %ld. Time per-call %g)\n\n",
      nsurf[dir], 1.0*tmsurf[dir]/CLOCKS_PER_SEC/nsurf[dir]);

    totSurfTm += 1.0*tmsurf[dir]/CLOCKS_PER_SEC;
  }

  printf("Total volume updates took %g. Surface updates took %g\n",
    totVolTm, totSurfTm);

  // release resources
  gkyl_array_release(fIn);
  gkyl_array_release(fOut);
  gkyl_array_release(em);
  gkyl_dg_eqn_release(vlasov_eqn);
  
  return 0;
}
