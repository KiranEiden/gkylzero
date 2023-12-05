#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gkgeom.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

#include <math.h>

static inline double sq(double x) { return x*x; }
static inline double cub(double x) { return x*x*x; }
static inline double qad(double x) { return x*x*x*x; }
static inline double pen(double x) { return x*x*x*x*x; }
static inline double hex(double x) { return x*x*x*x*x*x; }

static void
calc_exact_quad_2d(double dx, double dy, double xc, double yc, double coeff[16])
{
/* Exact expansion computed using the following Maxima code:

   load("basis-precalc/basisTensor2x")$
   load("modal-basis")$
   bc : basisC[3]$

   f : (1-x1^2)*y1^2$
   fS : subst([x1=xc+dx/2*x, y1=yc+dy/2*y],f)$
   proj : calcInnerProdList(varsC,1,bc,fS)$
 */
  
  coeff[0] = (-2.0*sq(xc)*sq(yc))-0.1666666666666667*sq(dx)*sq(yc)+2.0*sq(yc)-0.1666666666666667*sq(dy)*sq(xc)-0.01388888888888889*sq(dx)*sq(dy)+0.1666666666666667*sq(dy); 
  coeff[1] = (-1.154700538379252*dx*xc*sq(yc))-0.09622504486493764*dx*sq(dy)*xc; 
  coeff[2] = (-1.154700538379252*dy*sq(xc)*yc)-0.09622504486493764*sq(dx)*dy*yc+1.154700538379252*dy*yc; 
  coeff[3] = -0.6666666666666666*dx*dy*xc*yc; 
  coeff[4] = (-0.149071198499986*sq(dx)*sq(yc))-0.01242259987499883*sq(dx)*sq(dy); 
  coeff[5] = (-0.149071198499986*sq(dy)*sq(xc))-0.01242259987499883*sq(dx)*sq(dy)+0.149071198499986*sq(dy); 
  coeff[6] = -0.08606629658238703*sq(dx)*dy*yc; 
  coeff[7] = -0.08606629658238703*dx*sq(dy)*xc; 
  coeff[8] = 0.0; 
  coeff[9] = 0.0; 
  coeff[10] = -0.01111111111111111*sq(dx)*sq(dy); 
  coeff[11] = 0.0; 
  coeff[12] = 0.0; 
  coeff[13] = 0.0; 
  coeff[14] = 0.0; 
  coeff[15] = 0.0;
}

static void
cubic_1d(void)
{
  double lower[] = { 0.0 }, upper[] = { 10.0 };
  int cells[] = { 16 };

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  // nodal grid used in IO so we can plot things
  double nc_lower[] = { lower[0] - 0.5*grid.dx[0] };
  double nc_upper[] = { upper[0] + 0.5*grid.dx[0] };
  int nc_cells[] = { cells[0] + 1 };
  struct gkyl_rect_grid nc_grid;
  gkyl_rect_grid_init(&nc_grid, 1, nc_lower, nc_upper, nc_cells);

  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 0, 0 };  
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  struct gkyl_range nc_local, nc_local_ext;
  gkyl_create_grid_ranges(&nc_grid, nghost, &nc_local_ext, &nc_local);

  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, 1, 3);

  struct gkyl_array *psi_nodal = gkyl_array_new(GKYL_DOUBLE, 1, cells[0]+1);
  struct gkyl_array *psi_cubic = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_dg_basis_op_mem *mem = gkyl_dg_alloc_cubic_1d(cells[0]);
  double xn[1];
   
  do {
    // initialize 1D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);
      
      gkyl_rect_grid_ll_node(&grid, iter.idx, xn);
      
      double *pn = gkyl_array_fetch(psi_nodal, nidx);
      pn[0] = -sq(xn[0]) + 0.15*cub(xn[0]);
    }
    
    // compute cubic expansion
    gkyl_dg_calc_cubic_1d_from_nodal_vals(mem, cells[0], grid.dx[0],
      psi_nodal, psi_cubic);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, psi_nodal, "nodal_1d_a.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, psi_cubic, "cubic_1d_a.gkyl");
  } while (0);

  do {
    // initialize 1D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);
      
      gkyl_rect_grid_ll_node(&grid, iter.idx, xn);

      double *pn = gkyl_array_fetch(psi_nodal, nidx);
      pn[0] = sq(xn[0]) + 10*sq(sin(xn[0]));
    }
    // compute cubic expansion
    gkyl_dg_calc_cubic_1d_from_nodal_vals(mem, cells[0], grid.dx[0],
      psi_nodal, psi_cubic);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, psi_nodal, "nodal_1d_b.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, psi_cubic, "cubic_1d_b.gkyl");
  } while (0);  

  gkyl_array_release(psi_nodal);
  gkyl_array_release(psi_cubic);
  gkyl_dg_basis_op_mem_release(mem);
}

static void
cubic_2d(void)
{
  double lower[] = { 0.0, 0.0 }, upper[] = { 1.0, 1.0 };
  int cells[] = { 16, 8 };

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  // nodal grid used in IO so we can plot things
  double nc_lower[] = { lower[0] - 0.5*grid.dx[0], lower[1] - 0.5*grid.dx[1] };
  double nc_upper[] = { upper[0] + 0.5*grid.dx[0], upper[1] + 0.5*grid.dx[1] };
  int nc_cells[] = { cells[0] + 1, cells[1] + 1 };
  struct gkyl_rect_grid nc_grid;
  gkyl_rect_grid_init(&nc_grid, 2, nc_lower, nc_upper, nc_cells);

  struct gkyl_range local, local_ext;
  int nghost[GKYL_MAX_CDIM] = { 0, 0 };  
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  struct gkyl_range nc_local, nc_local_ext;
  gkyl_create_grid_ranges(&nc_grid, nghost, &nc_local_ext, &nc_local);

  struct gkyl_basis basis;
  gkyl_cart_modal_tensor(&basis, 2, 3);

  struct gkyl_array *psi_nodal = gkyl_array_new(GKYL_DOUBLE, 1, (cells[0]+1)*(cells[1]+1));
  struct gkyl_array *psi_cubic = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_dg_basis_op_mem *mem = gkyl_dg_alloc_cubic_2d(cells);

  double xn[2];
  
  do {
    // initialize 2D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);
      
      gkyl_rect_grid_ll_node(&grid, iter.idx, xn);
      
      double *pn = gkyl_array_fetch(psi_nodal, nidx);      
      pn[0] = (-sq(xn[0]) + 0.15*cub(xn[0]));
    }
    // compute cubic expansion
    gkyl_dg_calc_cubic_2d_from_nodal_vals(mem, cells, grid.dx,
      psi_nodal, psi_cubic);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, psi_nodal, "nodal_2d_a.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, psi_cubic, "cubic_2d_a.gkyl");
  } while (0);

  do {
    // initialize 2D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);

      gkyl_rect_grid_ll_node(&grid, iter.idx, xn);
      
      double *pn = gkyl_array_fetch(psi_nodal, nidx);      
      pn[0] = (-sq(xn[1]) + 0.15*cub(xn[1]));
    }
    // compute cubic expansion
    gkyl_dg_calc_cubic_2d_from_nodal_vals(mem, cells, grid.dx,
      psi_nodal, psi_cubic);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, psi_nodal, "nodal_2d_b.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, psi_cubic, "cubic_2d_b.gkyl");
  } while (0);

  do {
    // initialize 2D nodal values
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &nc_local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&nc_local, iter.idx);

      gkyl_rect_grid_ll_node(&grid, iter.idx, xn);

      double *pn = gkyl_array_fetch(psi_nodal, nidx);      
      pn[0] = (1-sq(xn[0]))*sq(xn[1]);
    }
    // compute cubic expansion
    gkyl_dg_calc_cubic_2d_from_nodal_vals(mem, cells, grid.dx,
      psi_nodal, psi_cubic);

    int ilo = local.lower[0], iup = local.upper[0];
    int jlo = local.lower[1], jup = local.upper[1];

    double err = 0.0;

    // check against exact expansion
    gkyl_range_iter_init(&iter, &local);
    while (gkyl_range_iter_next(&iter)) {
      long nidx = gkyl_range_idx(&local, iter.idx);
      
      double xc[2], coeff[16];
      gkyl_rect_grid_cell_center(&grid, iter.idx, xc);

      if ((iter.idx[0] != ilo) && (iter.idx[0] != iup) && (iter.idx[1] != jlo) && (iter.idx[1] != jup)) {
        calc_exact_quad_2d(grid.dx[0], grid.dx[1], xc[0], xc[1], coeff);
        
        const double *pc = gkyl_array_cfetch(psi_cubic, nidx);
        for (int i=0; i<16; ++i)
          err = fmax(err, fabs(pc[i]-coeff[i]));
      }
    }
    printf("Max error for a quadratic for interior cells is %lg\n", err);
    
    gkyl_grid_sub_array_write(&nc_grid, &nc_local, psi_nodal, "nodal_2d_c.gkyl");
    gkyl_grid_sub_array_write(&grid, &local, psi_cubic, "cubic_2d_c.gkyl");
  } while (0);

  gkyl_array_release(psi_nodal);
  gkyl_array_release(psi_cubic);
  gkyl_dg_basis_op_mem_release(mem);
}

int
main(int argc, char **argv)
{
  cubic_1d();
  cubic_2d();
  
  return 0;
}
    
