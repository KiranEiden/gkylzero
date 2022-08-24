#include <float.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_prop.h>

static inline double
fmax3(double a, double b, double c)
{
  return fmax(fmax(a,b),c);
}

static inline double
fmin3(double a, double b, double c)
{
  return fmin(fmin(a,b),c);
}

// limiter function
static inline double
limiter_function(double r, enum gkyl_wave_limiter limiter)
{
  double theta = 0.0;
  switch (limiter) {
    case GKYL_NO_LIMITER:
      theta = 1.0;
      break;
    
    case GKYL_MIN_MOD:
      theta = fmax(0, fmin(1, r));
      break;

    case GKYL_SUPERBEE:
      theta = fmax3(0.0, fmin(1, 2*r), fmin(2.0, r));
      break;

    case GKYL_VAN_LEER:
      theta = (r+fabs(r))/(1+fabs(r));
      break;

    case GKYL_MONOTONIZED_CENTERED:
      theta = fmax(0.0, fmin3((1.0+r)/2, 2, 2*r));
      break;

    case GKYL_BEAM_WARMING:
      theta = r;
      break;

    case GKYL_ZERO:
      theta = 0;
      break;
  }
  return theta;
}

gkyl_wave_prop*
gkyl_wave_prop_new(struct gkyl_wave_prop_inp winp)
{
  gkyl_wave_prop *up = gkyl_malloc(sizeof(gkyl_wave_prop));

  up->grid = *(winp.grid);
  up->ndim = up->grid.ndim;
  
  up->num_up_dirs = winp.num_up_dirs;
  for (int i=0; i<winp.num_up_dirs; ++i)
    up->update_dirs[i] = winp.update_dirs[i];

  up->limiter = winp.limiter == 0 ? GKYL_MONOTONIZED_CENTERED : winp.limiter;
  up->cfl = winp.cfl;
  up->equation = gkyl_wv_eqn_acquire(winp.equation);

  up->force_low_order_flux = winp.force_low_order_flux;
  up->check_inv_domain = winp.check_inv_domain;

  int nghost[3] = { 2, 2, 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&up->grid, nghost, &ext_range, &range);

  int max_1d = 0;
  for (int d=0; d<ext_range.ndim; ++d) {
    int shape = gkyl_range_shape(&ext_range, d);
    max_1d = max_1d > shape ? max_1d : shape;
  }

  // allocate memory to store 1D slices of waves, speeds and
  // second-order correction flux
  int meqn = winp.equation->num_equations, mwaves = winp.equation->num_waves;
  up->waves = gkyl_array_new(GKYL_DOUBLE, meqn*mwaves, max_1d);
  up->apdq = gkyl_array_new(GKYL_DOUBLE, meqn, max_1d);
  up->amdq = gkyl_array_new(GKYL_DOUBLE, meqn, max_1d);
  up->speeds = gkyl_array_new(GKYL_DOUBLE, mwaves, max_1d);
  up->flux2 = gkyl_array_new(GKYL_DOUBLE, meqn, max_1d);

  up->redo_fluct = gkyl_array_new(GKYL_DOUBLE, meqn, max_1d);

  // construct geometry
  up->geom = gkyl_wave_geom_acquire(winp.geom);

  up->n_calls = up->n_bad_advance_calls = 0;
  up->n_bad_cells = up->n_max_bad_cells = 0;

  return up;
}

// some helper functions

static inline void
copy_wv_vec(int n, double * GKYL_RESTRICT out, const double * GKYL_RESTRICT inp)
{
  for (int i=0; i<n; ++i) out[i] = inp[i];
}

static inline void
calc_jump(int n, const double *ql, const double *qr, double * GKYL_RESTRICT jump)
{
  for (int d=0; d<n; ++d) jump[d] = qr[d]-ql[d];
}

static inline void
calc_first_order_update(int meqn, double dtdx,
  double * GKYL_RESTRICT q, const double * GKYL_RESTRICT amdq_r, const double * GKYL_RESTRICT apdq_l)
{
  for (int i=0; i<meqn; ++i)
    q[i] = q[i] - dtdx*(apdq_l[i] + amdq_r[i]);
}

static inline double
calc_cfla(int mwaves, double cfla, double dtdx, const double *s)
{
  double c = cfla;
  for (int i=0; i<mwaves; ++i)
    c = fmax(c, dtdx*fabs(s[i]));
  return c;
}

static inline double
wave_dot_prod(int meqn, const double * GKYL_RESTRICT wa, const double * GKYL_RESTRICT wb)
{
  double dot = 0.0;
  for (int i=0; i<meqn; ++i) dot += wa[i]*wb[i];
  return dot;
}

static inline void
wave_rescale(int meqn, double fact, double *w)
{
  for (int i=0; i<meqn; ++i) w[i] *= fact; 
}

static inline void
calc_second_order_flux(int meqn, double dtdx, double s,
  const double *waves, double * GKYL_RESTRICT flux2)
{
  double sfact = 0.5*fabs(s)*(1-fabs(s)*dtdx);
  for (int i=0; i<meqn; ++i)
    flux2[i] += sfact*waves[i];
}

static inline void
calc_second_order_update(int meqn, double dtdx, double * GKYL_RESTRICT qout,
  const double *fl, const double *fr)
{
  for (int i=0; i<meqn; ++i)
    qout[i] += -dtdx*(fr[i]-fl[i]);
}

static void
limit_waves(const gkyl_wave_prop *wv, int mwaves,
  const struct gkyl_range *slice_range,
  int lower, int upper, struct gkyl_array *waves, const struct gkyl_array *speed)
{
  int meqn = wv->equation->num_equations;

  for (int mw=0; mw<mwaves; ++mw) {
    const double *wl = gkyl_array_cfetch(waves, gkyl_ridx(*slice_range, lower-1));
    const double *wr = gkyl_array_cfetch(waves, gkyl_ridx(*slice_range, lower));

    double dotr = wave_dot_prod(meqn, &wl[mw*meqn], &wr[mw*meqn]);

    for (int i=lower; i<=upper; ++i) {
      double dotl = dotr;
      
      double * GKYL_RESTRICT wi = gkyl_array_fetch(waves, gkyl_ridx(*slice_range, i));
      const double * GKYL_RESTRICT wi1 = gkyl_array_cfetch(waves, gkyl_ridx(*slice_range, i+1));
      
      double wnorm2 = wave_dot_prod(meqn, &wi[mw*meqn], &wi[mw*meqn]);
      dotr = wave_dot_prod(meqn, &wi[mw*meqn], &wi1[mw*meqn]);

      if (wnorm2 > 0) {
        const double *s = gkyl_array_cfetch(speed, gkyl_ridx(*slice_range, i));
        double r = s[mw] > 0 ? dotl/wnorm2 : dotr/wnorm2;
        double theta = limiter_function(r, wv->limiter);
        wave_rescale(meqn, theta, &wi[mw*meqn]);
      }
    }
  }
}

// advance method
struct gkyl_wave_prop_status
gkyl_wave_prop_advance(gkyl_wave_prop *wv,
  double tm, double dt, const struct gkyl_range *update_range,
  const struct gkyl_array *qin, struct gkyl_array *qout)
{
  wv->n_calls += 1;
  
  int ndim = update_range->ndim;
  int meqn = wv->equation->num_equations;
  //  when forced to use Lax fluxes, we only have a single wave
  int mwaves = wv->force_low_order_flux ? 1 :  wv->equation->num_waves;

  double cfla = 0.0, cfl = wv->cfl, cflm = 1.1*cfl;

  double ql_local[meqn], qr_local[meqn];
  double waves_local[meqn*mwaves];
  double amdq_local[meqn], apdq_local[meqn];
  double delta[meqn];

  int idxl[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];

  // state of the update
  enum update_state {
    WV_FIRST_SWEEP, WV_POSITIVITY_SWEEP, WV_FIN_SWEEP
  } state, next_state;

  for (int d=0; d<wv->num_up_dirs; ++d) {
    int dir = wv->update_dirs[d];

    double dtdx = dt/wv->grid.dx[dir];

    // upper/lower bounds in direction 'd'. These are edge indices
    int loidx = update_range->lower[dir]-1;
    int upidx = update_range->upper[dir]+2;

    // cell indices in 1D slice for interior cells
    int loidx_c = update_range->lower[dir];
    int upidx_c = update_range->upper[dir];

    struct gkyl_range slice_range;
    gkyl_range_init(&slice_range, 1, (int[]) { loidx }, (int[]) { upidx } );

    struct gkyl_range perp_range;
    gkyl_range_shorten(&perp_range, update_range, dir, 1);
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &perp_range);

    // outer loop is over perpendicular directions, inner loop over 1D
    // slice along that direction
    while (gkyl_range_iter_next(&iter)) {
      
      gkyl_copy_int_arr(ndim, iter.idx, idxl);
      gkyl_copy_int_arr(ndim, iter.idx, idxr);

      gkyl_array_clear(wv->redo_fluct, 1.0);
      
      enum gkyl_wv_flux_type ftype = wv->force_low_order_flux ?
        GKYL_WV_LOW_ORDER_FLUX : GKYL_WV_HIGH_ORDER_FLUX;

      state = WV_FIRST_SWEEP;

      // perform 1D sweeps, fixing positivity if required
      while (state != WV_FIN_SWEEP) {

        if (state == WV_POSITIVITY_SWEEP)
          ftype = GKYL_WV_LOW_ORDER_FLUX;

        // copy previous time-step solution
        for (int i=loidx_c; i<=upidx_c; ++i) {
          idxl[dir] = i; // cell index
          long lidx = gkyl_range_idx(update_range, idxl);
          copy_wv_vec(meqn, gkyl_array_fetch(qout, lidx), gkyl_array_cfetch(qin, lidx));
        }

        for (int i=loidx; i<upidx; ++i) { // CAUTION: should this be i<= as our indexing is inclusive?!
          idxl[dir] = i-1; idxr[dir] = i;
          long sidx = gkyl_ridx(slice_range, i);

          const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wv->geom, idxr);
          double *s = gkyl_array_fetch(wv->speeds, sidx);
          const double *redo_fluct = gkyl_array_cfetch(wv->redo_fluct, sidx);

          if (redo_fluct[0] > 0.0) {

            // compute fluctuations and waves only if needed (this
            // prevents doing the full 1D sweep with low-order fluxes
            // on positivity violations)
            long lidx = gkyl_range_idx(update_range, idxl);
            long ridx = gkyl_range_idx(update_range, idxr);

            const double *qinl = gkyl_array_cfetch(qin, lidx);
            const double *qinr = gkyl_array_cfetch(qin, ridx);

            wv->equation->rotate_to_local_func(cg->tau1[dir], cg->tau2[dir], cg->norm[dir], qinl, ql_local);
            wv->equation->rotate_to_local_func(cg->tau1[dir], cg->tau2[dir], cg->norm[dir], qinr, qr_local);

            calc_jump(meqn, ql_local, qr_local, delta);
          

            gkyl_wv_eqn_waves(wv->equation, ftype, delta, ql_local, qr_local, waves_local, s);

            double lenr = cg->lenr[dir];
            for (int mw=0; mw<mwaves; ++mw)
              s[mw] *= lenr; // rescale speeds

            // compute fluctuations in local coordinates
            gkyl_wv_eqn_qfluct(wv->equation, ftype, ql_local, qr_local,
              waves_local, s, amdq_local, apdq_local);
        
            double *waves = gkyl_array_fetch(wv->waves, sidx);
            for (int mw=0; mw<mwaves; ++mw)
              // rotate waves back
              wv->equation->rotate_to_global_func(
                cg->tau1[dir], cg->tau2[dir], cg->norm[dir], &waves_local[mw*meqn], &waves[mw*meqn]
              );

            // rotate fluctuations
            double *amdq = gkyl_array_fetch(wv->amdq, sidx);
            wv->equation->rotate_to_global_func(
              cg->tau1[dir], cg->tau2[dir], cg->norm[dir], amdq_local, amdq);
            
            double *apdq = gkyl_array_fetch(wv->apdq, sidx);
            wv->equation->rotate_to_global_func(
              cg->tau1[dir], cg->tau2[dir], cg->norm[dir], apdq_local, apdq);
          }
          
          cfla = calc_cfla(mwaves, cfla, dtdx/cg->kappa, s);
        }

        if (cfla > cflm) // check time-step before any updates are performed
          return (struct gkyl_wave_prop_status) { .success = 0, .dt_suggested = dt*cfl/cfla };

        // compute first-order update in each cell
        for (int i=loidx_c; i<=upidx_c; ++i) { // loop is over cells
        
          idxl[dir] = i; // cell index and left-edge index
          long lidx = gkyl_range_idx(update_range, idxl);

          const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wv->geom, idxl);

          calc_first_order_update(meqn, dtdx/cg->kappa,
            gkyl_array_fetch(qout, lidx), 
            gkyl_array_cfetch(wv->amdq, gkyl_ridx(slice_range, i+1)),
            gkyl_array_cfetch(wv->apdq, gkyl_ridx(slice_range, i))
          );
        }

        if (state == WV_FIRST_SWEEP) {
          // we only compute second-correction if we are in first sweep
          
          // apply limiters to waves for all edges in update range,
          // including edges that are on the range boundary
          limit_waves(wv, mwaves, &slice_range,
            update_range->lower[dir], update_range->upper[dir]+1, wv->waves, wv->speeds);

          // get the kappa in the first ghost cell on left (needed in
          // the second order flux calculation)
          idxl[dir] = update_range->lower[dir]-1;
          const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wv->geom, idxl);
          double kappal = cg->kappa;

          gkyl_array_clear(wv->flux2, 0.0);
          // compute second-order correction fluxes at each interface:
          // note that there is one extra edge than cell
          for (int i=loidx_c; i<=upidx_c+1; ++i) {
            long sidx = gkyl_ridx(slice_range, i);

            const double *waves = gkyl_array_cfetch(wv->waves, sidx);
            const double *s = gkyl_array_cfetch(wv->speeds, sidx);
            double *flux2 = gkyl_array_fetch(wv->flux2, sidx);

            idxl[dir] = i;
            const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wv->geom, idxl);
            double kappar = cg->kappa;

            for (int mw=0; mw<mwaves; ++mw)
              calc_second_order_flux(meqn, dtdx/(0.5*(kappal+kappar)), s[mw], &waves[mw*meqn], flux2);

            kappal = kappar;
          }

          // add second correction flux to solution in each interior cell
          for (int i=loidx_c; i<=upidx_c; ++i) {
            long sidx = gkyl_ridx(slice_range, i);

            idxl[dir] = i;
            const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wv->geom, idxl);

            calc_second_order_update(meqn, dtdx/cg->kappa,
              gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl)),
              gkyl_array_cfetch(wv->flux2, gkyl_ridx(slice_range, i)),
              gkyl_array_cfetch(wv->flux2, gkyl_ridx(slice_range, i+1))
            );
          }
        }

        next_state = WV_FIN_SWEEP;
        // check invariant domains if needed
        if ( (state == WV_FIRST_SWEEP) && wv->check_inv_domain) {
          long n_bad_cells = 0;            

          gkyl_array_clear(wv->redo_fluct, 0.0); // by default no edge needs recomputing
          
          // check if invariant domains are violated, flagging edges
          // of each bad cell
          for (int i=loidx_c; i<=upidx_c; ++i) {
            idxl[dir] = i;
            const double *qt = gkyl_array_cfetch(qout, gkyl_range_idx(update_range, idxl));
            if (!wv->equation->check_inv_func(wv->equation, qt)) {

              double *redo_fluct_l = gkyl_array_fetch(wv->redo_fluct, gkyl_ridx(slice_range, i));
              double *redo_fluct_r = gkyl_array_fetch(wv->redo_fluct, gkyl_ridx(slice_range, i+1));
              // mark left and right edges so fluctuations are redone
              redo_fluct_l[0] = 1.0;
              redo_fluct_r[0] = 1.0;

              n_bad_cells += 1;
            }
          }

          if (n_bad_cells > 0) {
            // we need to resweep the 1D slice again
            next_state = WV_POSITIVITY_SWEEP;
            wv->n_bad_advance_calls += 1;
          }

          wv->n_bad_cells += n_bad_cells;
          wv->n_max_bad_cells = wv->n_max_bad_cells >  n_bad_cells ? wv->n_max_bad_cells : n_bad_cells;
        }

        state = next_state; // change state for next sweep
        
      } // end loop over sweeps
    } // end loop over perpendicular directions
  } // end loop over directions

  // compute allowable time-step from this update, but suggest only
  // bigger time-step; (Only way dt can reduce is if the update
  // fails. If the code comes here the update suceeded and so we
  // should not allow dt to reduce).
  double dt_suggested = dt*cfl/fmax(cfla, DBL_MIN);

  return (struct gkyl_wave_prop_status) {
    .success = 1,
    .dt_suggested = dt_suggested > dt ? dt_suggested : dt
  };
}

double
gkyl_wave_prop_max_dt(const gkyl_wave_prop *wv, const struct gkyl_range *update_range,
  const struct gkyl_array *qin)
{
  double max_dt = DBL_MAX;
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {

    for (int d=0; d<wv->num_up_dirs; ++d) {
      int dir = wv->update_dirs[d];
      double dx = wv->grid.dx[dir];

      const double *q = gkyl_array_cfetch(qin, gkyl_range_idx(update_range, iter.idx));
      double maxs = gkyl_wv_eqn_max_speed(wv->equation, q);
      max_dt = fmin(max_dt, wv->cfl*dx/maxs);
    }
    
  }

  return max_dt;
}

struct gkyl_wave_prop_stats
gkyl_wave_prop_stats(const gkyl_wave_prop *wv)
{
  return (struct gkyl_wave_prop_stats) {
    .n_calls = wv->n_calls,
    .n_bad_advance_calls = wv->n_bad_advance_calls,
    .n_bad_cells = wv->n_bad_cells,
    .n_max_bad_cells = wv->n_max_bad_cells
  };
}

void
gkyl_wave_prop_release(gkyl_wave_prop* up)
{
  gkyl_wv_eqn_release(up->equation);
  gkyl_array_release(up->waves);
  gkyl_array_release(up->apdq);
  gkyl_array_release(up->amdq);
  gkyl_array_release(up->speeds);
  gkyl_array_release(up->flux2);
  gkyl_array_release(up->redo_fluct);
  
  gkyl_wave_geom_release(up->geom);
  
  gkyl_free(up);
}
