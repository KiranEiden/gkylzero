#include <gkyl_fem_poisson.h>
#include <gkyl_fem_poisson_priv.h>

gkyl_fem_poisson*
gkyl_fem_poisson_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis basis,
  struct gkyl_poisson_bc *bcs, double epsilon_const, struct gkyl_array *epsilon_var, struct gkyl_array *kSq, bool use_gpu)
{

  gkyl_fem_poisson *up = gkyl_malloc(sizeof(gkyl_fem_poisson));

  up->kernels = gkyl_malloc(sizeof(struct gkyl_fem_poisson_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_fem_poisson_kernels));
  else
    up->kernels_cu = up->kernels;
#else
  up->kernels_cu = up->kernels;
#endif

  up->ndim = grid->ndim;
  up->grid = *grid;
  up->num_basis =  basis.num_basis;
  up->basis_type = basis.b_type;
  up->poly_order = basis.poly_order;
  up->basis = basis;
  up->use_gpu = use_gpu;

  if (epsilon_var) {
    up->isvareps = true;
    up->epsilon  = epsilon_var;
  } else {
    up->isvareps = false;
    // Create an array to hold epsilon.
#ifdef GKYL_HAVE_CUDA
    up->epsilon = up->use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 1) : gkyl_array_new(GKYL_DOUBLE, 1, 1);
#else
    up->epsilon = gkyl_array_new(GKYL_DOUBLE, 1, 1);
#endif
    gkyl_array_shiftc(up->epsilon, epsilon_const, 0);
  }

  // We assume epsilon_var and kSq live on the device, and we create a host-side
  // copies temporarily to compute the LHS matrix. This also works for CPU solves.
  struct gkyl_array *epsilon_ho = gkyl_array_new(GKYL_DOUBLE, up->epsilon->ncomp, up->epsilon->size);
  gkyl_array_copy(epsilon_ho, up->epsilon);
  struct gkyl_array *kSq_ho;
  if (kSq) {
    up->ishelmholtz = true;
    kSq_ho = gkyl_array_new(GKYL_DOUBLE, kSq->ncomp, kSq->size);
    gkyl_array_copy(kSq_ho, kSq);
  } else {
    up->ishelmholtz = false;
    kSq_ho = gkyl_array_new(GKYL_DOUBLE, 1, 1);
    gkyl_array_clear(kSq_ho, 0.);
  }

  up->globalidx = gkyl_malloc(sizeof(long[up->num_basis])); // global index, one for each basis in a cell.

  // Local and local-ext ranges for whole-grid arrays.
  int ghost[GKYL_MAX_DIM];
  for (int d=0; d<up->ndim; d++) ghost[d] = 1;
  gkyl_create_grid_ranges(grid, ghost, &up->local_range_ext, &up->local_range);
  // Range of cells we'll solve Poisson in, as
  // a sub-range of up->local_range_ext.
  int sublower[POISSON_MAX_DIM], subupper[POISSON_MAX_DIM];
  for (int d=0; d<up->ndim; d++) {
    sublower[d] = up->local_range.lower[d];
    subupper[d] = up->local_range.upper[d];
  }
  gkyl_sub_range_init(&up->solve_range, &up->local_range_ext, sublower, subupper);
  for (int d=0; d<up->ndim; d++) up->num_cells[d] = up->solve_range.upper[d]-up->solve_range.lower[d]+1;

  // Prepare for periodic domain case.
  for (int d=0; d<up->ndim; d++) {
    // Sanity check.
    if ((bcs->lo_type[d] == GKYL_POISSON_PERIODIC && bcs->up_type[d] != GKYL_POISSON_PERIODIC) ||
        (bcs->lo_type[d] != GKYL_POISSON_PERIODIC && bcs->up_type[d] == GKYL_POISSON_PERIODIC))
      assert(false);
  }
  for (int d=0; d<up->ndim; d++) up->isdirperiodic[d] = bcs->lo_type[d] == GKYL_POISSON_PERIODIC;
  up->isdomperiodic = true;
  for (int d=0; d<up->ndim; d++) up->isdomperiodic = up->isdomperiodic && up->isdirperiodic[d];
  if (up->isdomperiodic) {
#ifdef GKYL_HAVE_CUDA
    if (up->use_gpu) {
      up->rhs_cellavg = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, up->local_range_ext.volume);
      up->rhs_avg_cu = (double*) gkyl_cu_malloc(sizeof(double)); 
    } else {
      up->rhs_cellavg = gkyl_array_new(GKYL_DOUBLE, 1, up->local_range_ext.volume);
    }
#else
    up->rhs_cellavg = gkyl_array_new(GKYL_DOUBLE, 1, up->local_range_ext.volume);
#endif
    up->rhs_avg = (double*) gkyl_malloc(sizeof(double)); 
    gkyl_array_clear(up->rhs_cellavg, 0.0);
    // Factor accounting for normalization when subtracting a constant from a
    // DG field and the 1/N to properly compute the volume averaged RHS.
    up->mavgfac = -pow(sqrt(2.),up->ndim)/up->solve_range.volume;
  }

  // Pack BC values into a single array for easier use in kernels.
  for (int d=0; d<up->ndim; d++) {
    for (int k=0; k<6; k++) up->bcvals[d*2*3+k] = 0.0; // default. Not used in some cases (e.g. periodic).
    if (bcs->lo_type[d] != GKYL_POISSON_PERIODIC) {
      int vnum, voff;
      vnum = bcs->lo_type[d] == GKYL_POISSON_ROBIN ? 3 : 1;
      voff = bcs->lo_type[d] == GKYL_POISSON_ROBIN ? 0 : 2;
      for (int k=0; k<vnum; k++) up->bcvals[d*2*3+voff+k] = bcs->lo_value[d].v[k];

      vnum = bcs->up_type[d] == GKYL_POISSON_ROBIN ? 3 : 1;
      voff = bcs->up_type[d] == GKYL_POISSON_ROBIN ? 0 : 2;
      for (int k=0; k<vnum; k++) up->bcvals[d*2*3+voff+3+k] = bcs->up_value[d].v[k];
    }
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    up->bcvals_cu = (double *) gkyl_cu_malloc(sizeof(double[POISSON_MAX_DIM*3*2]));
    gkyl_cu_memcpy(up->bcvals_cu, up->bcvals, sizeof(double[POISSON_MAX_DIM*3*2]), GKYL_CU_MEMCPY_H2D);
  }
#endif

  // Compute the number of local and global nodes.
  up->numnodes_local = up->num_basis;
  up->numnodes_global = gkyl_fem_poisson_global_num_nodes(up->ndim, up->poly_order, basis.b_type, up->num_cells, up->isdirperiodic);

  for (int d=0; d<up->ndim; d++) up->dx[d] = up->grid.dx[d];  // Cell lengths.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    up->dx_cu = (double *) gkyl_cu_malloc(sizeof(double[POISSON_MAX_DIM]));
    gkyl_cu_memcpy(up->dx_cu, up->dx, sizeof(double[POISSON_MAX_DIM]), GKYL_CU_MEMCPY_H2D);
  }
#endif

  up->brhs = gkyl_array_new(GKYL_DOUBLE, 1, up->numnodes_global); // Global right side vector.

  // Select local-to-global mapping kernels:
  fem_poisson_choose_local2global_kernels(&basis, up->isdirperiodic, up->kernels->l2g);

  // Select lhs kernels:
  fem_poisson_choose_lhs_kernels(&basis, bcs, up->isvareps, up->kernels->lhsker);

  // Select rhs src kernels:
  fem_poisson_choose_src_kernels(&basis, bcs, up->isvareps, up->kernels->srcker);

  // Select sol kernel:
  up->kernels->solker = fem_poisson_choose_sol_kernels(&basis);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    fem_poisson_choose_kernels_cu(&basis, bcs, up->isvareps, up->isdirperiodic, up->kernels_cu);
#endif

  // Create a linear Ax=B problem. Here A is the discrete (global) matrix
  // representation of the LHS of the Helmholtz equation.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) 
    up->prob_cu = gkyl_cusolver_prob_new(up->numnodes_global, up->numnodes_global, 1);
  else
    up->prob = gkyl_superlu_prob_new(up->numnodes_global, up->numnodes_global, 1);
#else
  up->prob = gkyl_superlu_prob_new(up->numnodes_global, up->numnodes_global, 1);
#endif

  // Assign non-zero elements in A.
  gkyl_mat_triples *tri = gkyl_mat_triples_new(up->numnodes_global, up->numnodes_global);
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) gkyl_mat_triples_set_rowmaj_order(tri);
#endif
  gkyl_range_iter_init(&up->solve_iter, &up->solve_range);
  int idx0[POISSON_MAX_DIM];
  while (gkyl_range_iter_next(&up->solve_iter)) {
    long linidx = gkyl_range_idx(&up->solve_range, up->solve_iter.idx);

    double *eps_p = up->isvareps? gkyl_array_fetch(epsilon_ho, linidx) : epsilon_ho->data;
    double *kSq_p = up->ishelmholtz? gkyl_array_fetch(kSq_ho, linidx) : kSq_ho->data;

    int keri = idx_to_inup_ker(up->ndim, up->num_cells, up->solve_iter.idx);
    for (size_t d=0; d<up->ndim; d++) idx0[d] = up->solve_iter.idx[d]-1;
    up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

    // Apply the -nabla . (epsilon*nabla)-kSq stencil.
    keri = idx_to_inloup_ker(up->ndim, up->num_cells, up->solve_iter.idx);
    up->kernels->lhsker[keri](eps_p, kSq_p, up->dx, up->bcvals, up->globalidx, tri);
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cusolver_amat_from_triples(up->prob_cu, tri);
  } else {
    gkyl_superlu_amat_from_triples(up->prob, tri);
  }
#else
  gkyl_superlu_amat_from_triples(up->prob, tri);
#endif

  gkyl_mat_triples_release(tri);
  gkyl_array_release(epsilon_ho);
  gkyl_array_release(kSq_ho);

  return up;
}

void
gkyl_fem_poisson_set_rhs(gkyl_fem_poisson* up, struct gkyl_array *rhsin)
{

  if (up->isdomperiodic && !(up->ishelmholtz)) {
    // Subtract the volume averaged RHS from the RHS.
    gkyl_array_clear(up->rhs_cellavg, 0.0);
    gkyl_dg_calc_average_range(up->basis, 0, up->rhs_cellavg, 0, rhsin, up->solve_range);
#ifdef GKYL_HAVE_CUDA
    if (up->use_gpu) {
      gkyl_array_reduce_range(up->rhs_avg_cu, up->rhs_cellavg, GKYL_SUM, up->solve_range);
      gkyl_cu_memcpy(up->rhs_avg, up->rhs_avg_cu, sizeof(double), GKYL_CU_MEMCPY_D2H);
    } else {
      gkyl_array_reduce_range(up->rhs_avg, up->rhs_cellavg, GKYL_SUM, up->solve_range);
    }
#else
    gkyl_array_reduce_range(up->rhs_avg, up->rhs_cellavg, GKYL_SUM, up->solve_range);
#endif
    gkyl_array_shiftc(rhsin, up->mavgfac*up->rhs_avg[0], 0);
  }


#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(rhsin));

    gkyl_fem_poisson_set_rhs_cu(up, rhsin);
    return;
  }
#endif

  gkyl_array_clear(up->brhs, 0.0);

  gkyl_range_iter_init(&up->solve_iter, &up->solve_range);
  int idx0[POISSON_MAX_DIM];
  double *brhs_p = gkyl_array_fetch(up->brhs, 0);
  while (gkyl_range_iter_next(&up->solve_iter)) {
    long linidx = gkyl_range_idx(&up->solve_range, up->solve_iter.idx);

    double *eps_p = up->isvareps? gkyl_array_fetch(up->epsilon, linidx) : up->epsilon->data;
    double *rhsin_p = gkyl_array_fetch(rhsin, linidx);

    int keri = idx_to_inup_ker(up->ndim, up->num_cells, up->solve_iter.idx);
    for (size_t d=0; d<up->ndim; d++) idx0[d] = up->solve_iter.idx[d]-1;
    up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

    // Apply the RHS source stencil. It's mostly the mass matrix times a
    // modal-to-nodal operator times the source, modified by BCs in skin cells.
    keri = idx_to_inloup_ker(up->ndim, up->num_cells, up->solve_iter.idx);
    up->kernels->srcker[keri](eps_p, up->dx, rhsin_p, up->bcvals, up->globalidx, brhs_p);
  }

  gkyl_superlu_brhs_from_array(up->prob, brhs_p);

}

void
gkyl_fem_poisson_solve(gkyl_fem_poisson* up, struct gkyl_array *phiout) {
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(phiout));
    gkyl_fem_poisson_solve_cu(up, phiout);
    return;
  }
#endif

  gkyl_superlu_solve(up->prob);

  gkyl_range_iter_init(&up->solve_iter, &up->solve_range);
  int idx0[POISSON_MAX_DIM];
  gkyl_array_clear(phiout, 0.0);
  while (gkyl_range_iter_next(&up->solve_iter)) {
    long linidx = gkyl_range_idx(&up->solve_range, up->solve_iter.idx);

    double *phiout_p = gkyl_array_fetch(phiout, linidx);

    int keri = idx_to_inup_ker(up->ndim, up->num_cells, up->solve_iter.idx);
    for (size_t d=0; d<up->ndim; d++) idx0[d] = up->solve_iter.idx[d]-1;
    up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

    up->kernels->solker(gkyl_superlu_get_rhs_ptr(up->prob, 0), up->globalidx, phiout_p);

  }

}

void gkyl_fem_poisson_release(gkyl_fem_poisson *up)
{
  if (up->isdomperiodic) {
    gkyl_array_release(up->rhs_cellavg);
    gkyl_free(up->rhs_avg);
  }

  if (!up->isvareps)
    gkyl_array_release(up->epsilon);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->kernels_cu);
    gkyl_cu_free(up->dx_cu);
    if (up->isdomperiodic) gkyl_cu_free(up->rhs_avg_cu);
    gkyl_cu_free(up->bcvals_cu);
    gkyl_cusolver_prob_release(up->prob_cu);
  } else {
    gkyl_superlu_prob_release(up->prob);
  }
#else
  gkyl_superlu_prob_release(up->prob);
#endif

  gkyl_free(up->globalidx);
  gkyl_free(up->brhs);
  gkyl_free(up->kernels);
  gkyl_free(up);
}
