
#include <gkyl_vlasov_priv.h>

gkyl_vlasov_app*
gkyl_vlasov_app_new(struct gkyl_vm vm)
{
  assert(vm.num_species <= GKYL_MAX_SPECIES);
  
  gkyl_vlasov_app *app = gkyl_malloc(sizeof(gkyl_vlasov_app));
  
  int cdim = app->cdim = vm.cdim;
  int vdim = app->vdim = vm.vdim;
  int pdim = cdim+vdim;
  int poly_order = app->poly_order = vm.poly_order;
  int ns = app->num_species = vm.num_species;

  double cfl_frac = vm.cfl_frac == 0 ? 1.0 : vm.cfl_frac;
  app->cfl = cfl_frac/(2*poly_order+1);

#ifdef GKYL_HAVE_CUDA
  app->use_gpu = vm.use_gpu;
#else
  app->use_gpu = false; // can't use GPUs if we don't have them!
#endif
  
  app->num_periodic_dir = vm.num_periodic_dir;
  for (int d=0; d<cdim; ++d)
    app->periodic_dirs[d] = vm.periodic_dirs[d];

  strcpy(app->name, vm.name);
  app->tcurr = 0.0; // reset on init

  // check if there is a job pool
  if (vm.job_pool)
    app->job_pool = gkyl_job_pool_acquire(vm.job_pool);
  else
    app->job_pool = gkyl_null_pool_new(1);

  // basis functions
  switch (vm.basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      gkyl_cart_modal_serendip(&app->basis, pdim, poly_order);
      gkyl_cart_modal_serendip(&app->confBasis, cdim, poly_order);
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      gkyl_cart_modal_tensor(&app->basis, pdim, poly_order);
      gkyl_cart_modal_tensor(&app->confBasis, cdim, poly_order);
      break;

    default:
      assert(false);
      break;
  }

  gkyl_rect_grid_init(&app->grid, cdim, vm.lower, vm.upper, vm.cells);

  int ghost[] = { 1, 1, 1 };  
  gkyl_create_grid_ranges(&app->grid, ghost, &app->local_ext, &app->local);
  skin_ghost_ranges_init(&app->skin_ghost, &app->local_ext, ghost);

  app->has_field = !vm.skip_field; // note inversion of truth value
  
  if (app->has_field)
    app->field = vm_field_new(&vm, app);

  // allocate space to store species objects
  app->species = ns>0 ? gkyl_malloc(sizeof(struct vm_species[ns])) : 0;
  // create species grid & ranges
  for (int i=0; i<ns; ++i) {
    app->species[i].info = vm.species[i];
    vm_species_init(&vm, app, &app->species[i]);
  }

  // initialize stat object
  app->stat = (struct gkyl_vlasov_stat) {
    .use_gpu = app->use_gpu,
    .stage_2_dt_diff = { DBL_MAX, 0.0 },
    .stage_3_dt_diff = { DBL_MAX, 0.0 },
  };
  
  return app;
}

void
gkyl_vlasov_app_apply_ic(gkyl_vlasov_app* app, double t0)
{
  app->tcurr = t0;
  if (app->has_field)
    gkyl_vlasov_app_apply_ic_field(app, t0);
  for (int i=0;  i<app->num_species; ++i)
    gkyl_vlasov_app_apply_ic_species(app, i, t0);
}

void
gkyl_vlasov_app_apply_ic_field(gkyl_vlasov_app* app, double t0)
{
  app->tcurr = t0;

  struct timespec wtm = gkyl_wall_clock();
  vm_field_apply_ic(app, app->field, t0);
  app->stat.init_field_tm += gkyl_time_diff_now_sec(wtm);
    
  vm_field_apply_bc(app, app->field, app->field->em);  
}

void
gkyl_vlasov_app_apply_ic_species(gkyl_vlasov_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);

  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();
  vm_species_apply_ic(app, &app->species[sidx], t0);
  app->stat.init_species_tm += gkyl_time_diff_now_sec(wtm);

  vm_species_apply_bc(app, &app->species[sidx], app->species[sidx].f);
}

void
gkyl_vlasov_app_calc_mom(gkyl_vlasov_app* app)
{
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *s = &app->species[i];
    
    for (int m=0; m<app->species[i].info.num_diag_moments; ++m) {
      struct timespec wst = gkyl_wall_clock();
      vm_species_moment_calc(&s->moms[m], s->local, app->local, s->f);
      app->stat.mom_tm += gkyl_time_diff_now_sec(wst);
      app->stat.nmom += 1;
    }
  }
}

void
gkyl_vlasov_app_write(gkyl_vlasov_app* app, double tm, int frame)
{
  if (app->has_field)
    gkyl_vlasov_app_write_field(app, tm, frame);
  for (int i=0; i<app->num_species; ++i)
    gkyl_vlasov_app_write_species(app, i, tm, frame);
}
 
void
gkyl_vlasov_app_write_field(gkyl_vlasov_app* app, double tm, int frame)
{
  const char *fmt = "%s-field_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, frame);

  if (app->use_gpu) {
    // copy data from device to host before writing it out    
    gkyl_array_copy(app->field->em_host, app->field->em);
    gkyl_grid_sub_array_write(&app->grid, &app->local, app->field->em_host, fileNm);
  }
  else {
    gkyl_grid_sub_array_write(&app->grid, &app->local, app->field->em, fileNm);
  }
}

void
gkyl_vlasov_app_write_species(gkyl_vlasov_app* app, int sidx, double tm, int frame)
{
  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, app->species[sidx].info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[sidx].info.name, frame);

  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(app->species[sidx].f_host, app->species[sidx].f);
    gkyl_grid_sub_array_write(&app->species[sidx].grid, &app->species[sidx].local,
      app->species[sidx].f_host, fileNm);
  }
  else {
    gkyl_grid_sub_array_write(&app->species[sidx].grid, &app->species[sidx].local,
      app->species[sidx].f, fileNm);
  }
}

void
gkyl_vlasov_app_write_mom(gkyl_vlasov_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {
    for (int m=0; m<app->species[i].info.num_diag_moments; ++m) {
      
      const char *fmt = "%s-%s-%s_%d.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, app->species[i].info.name,
        app->species[i].info.diag_moments[m], frame);
      char fileNm[sz+1]; // ensures no buffer overflow  
      snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[i].info.name,
        app->species[i].info.diag_moments[m], frame);

      if (app->use_gpu)
        gkyl_array_copy(app->species[i].moms[m].marr_host, app->species[i].moms[m].marr);
      
      gkyl_grid_sub_array_write(&app->grid, &app->local, app->species[i].moms[m].marr_host, fileNm);
    }
  }
}

// Take a forward Euler step with the suggested time-step dt. This may
// not be the actual time-step taken. However, the function will never
// take a time-step larger than dt even if it is allowed by
// stability. The actual time-step and dt_suggested are returned in
// the status object.
static void
forward_euler(gkyl_vlasov_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], const struct gkyl_array *emin,
  struct gkyl_array *fout[], struct gkyl_array *emout, struct gkyl_update_status *st)
{
  app->stat.nfeuler += 1;
  
  double dtmin = DBL_MAX;

  // compute RHS of Vlasov equations
  for (int i=0; i<app->num_species; ++i) {
    if (app->has_field) {
      double qbym = app->species[i].info.charge/app->species[i].info.mass;
      gkyl_array_set(app->field->qmem, qbym, emin);
    }
    
    double dt1 = vm_species_rhs(app, &app->species[i], fin[i], app->has_field ? app->field->qmem : 0, fout[i]);
    dtmin = fmin(dtmin, dt1);
  }
  // compute RHS of Maxwell equations
  if (app->has_field) {
    double dt1 = vm_field_rhs(app, app->field, emin, emout);
    dtmin = fmin(dtmin, dt1);
  }

  double dt_max_rel_diff = 0.01;
  // check if dtmin is slightly smaller than dt. Use dt if it is
  // (avoids retaking steps if dt changes are very small).
  double dt_rel_diff = (dt-dtmin)/dt;
  if (dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff)
    dtmin = dt;

  // don't take a time-step larger that input dt
  double dta = st->dt_actual = dt < dtmin ? dt : dtmin;
  st->dt_suggested = dtmin;

  // complete update of distribution function
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_accumulate_range(gkyl_array_scale_range(fout[i], dta, app->species[i].local),
      1.0, fin[i], app->species[i].local);
    vm_species_apply_bc(app, &app->species[i], fout[i]);
  }

  if (app->has_field) {
    struct timespec wst = gkyl_wall_clock();
    // accumulate current contribution to electric field terms
    for (int i=0; i<app->num_species; ++i) {
      struct vm_species *s = &app->species[i];
      vm_species_moment_calc(&s->m1i, s->local, app->local, fin[i]);
    
      double qbyeps = s->info.charge/app->field->info.epsilon0;
      gkyl_array_accumulate_range(emout, -qbyeps, s->m1i.marr, app->local);
    }
    app->stat.current_tm += gkyl_time_diff_now_sec(wst);
  
    // complete update of field
    gkyl_array_accumulate_range(gkyl_array_scale_range(emout, dta, app->local),
      1.0, emin, app->local);
    vm_field_apply_bc(app, app->field, emout);
  }
}

// Take time-step using the RK3 method. Also sets the status object
// which has the actual and suggested dts used. These can be different
// from the actual time-step.
static struct gkyl_update_status
rk3(gkyl_vlasov_app* app, double dt0)
{
  const struct gkyl_array *fin[app->num_species];
  struct gkyl_array *fout[app->num_species];
  struct gkyl_update_status st = { .success = true };

  // time-stepper state
  enum { RK_STAGE_1, RK_STAGE_2, RK_STAGE_3, RK_COMPLETE } state = RK_STAGE_1;

  double tcurr = app->tcurr, dt = dt0;
  while (state != RK_COMPLETE) {
    switch (state) {
      case RK_STAGE_1:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f;
          fout[i] = app->species[i].f1;
        }
        forward_euler(app, tcurr, dt, fin, app->has_field ? app->field->em : 0,
          fout, app->has_field ? app->field->em1 : 0,
          &st
        );
        dt = st.dt_actual;
        state = RK_STAGE_2;
        break;

      case RK_STAGE_2:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f1;
          fout[i] = app->species[i].fnew;
        }
        forward_euler(app, tcurr+dt, dt, fin, app->has_field ? app->field->em1 : 0,
          fout, app->has_field ? app->field->emnew : 0,
          &st
        );
        if (st.dt_actual < dt) {

          // collect stats
          double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
          app->stat.stage_2_dt_diff[0] = fmin(app->stat.stage_2_dt_diff[0],
            dt_rel_diff);
          app->stat.stage_2_dt_diff[1] = fmax(app->stat.stage_2_dt_diff[1],
            dt_rel_diff);
          app->stat.nstage_2_fail += 1;
            
          dt = st.dt_actual;
          state = RK_STAGE_1; // restart from stage 1

        }
        else {
          for (int i=0; i<app->num_species; ++i)
            array_combine(app->species[i].f1,
              3.0/4.0, app->species[i].f, 1.0/4.0, app->species[i].fnew, app->species[i].local_ext);
          if (app->has_field)
            array_combine(app->field->em1,
              3.0/4.0, app->field->em, 1.0/4.0, app->field->emnew, app->local_ext);

          state = RK_STAGE_3;
        }
        break;

      case RK_STAGE_3:
        for (int i=0; i<app->num_species; ++i) {
          fin[i] = app->species[i].f1;
          fout[i] = app->species[i].fnew;
        }
        forward_euler(app, tcurr+dt/2, dt, fin, app->has_field ? app->field->em1 : 0,
          fout, app->has_field ? app->field->emnew : 0,
          &st
        );
        if (st.dt_actual < dt) {
          // collect stats
          double dt_rel_diff = (dt-st.dt_actual)/st.dt_actual;
          app->stat.stage_3_dt_diff[0] = fmin(app->stat.stage_3_dt_diff[0],
            dt_rel_diff);
          app->stat.stage_3_dt_diff[1] = fmax(app->stat.stage_3_dt_diff[1],
            dt_rel_diff);
          app->stat.nstage_3_fail += 1;
            
          dt = st.dt_actual;
          state = RK_STAGE_1; // restart from stage 1
            
          app->stat.nstage_2_fail += 1;
        }
        else {
          for (int i=0; i<app->num_species; ++i) {
            array_combine(app->species[i].f1,
              1.0/3.0, app->species[i].f, 2.0/3.0, app->species[i].fnew, app->species[i].local_ext);
            gkyl_array_copy_range(app->species[i].f, app->species[i].f1, app->species[i].local_ext);
          }
          if (app->has_field) {
            array_combine(app->field->em1,
              1.0/3.0, app->field->em, 2.0/3.0, app->field->emnew, app->local_ext);
            gkyl_array_copy_range(app->field->em, app->field->em1, app->local_ext);
          }

          state = RK_COMPLETE;
        }
        break;

      case RK_COMPLETE: // can't happen: suppresses warning
        break;
    }
  }
  
  return st;
}

struct gkyl_update_status
gkyl_vlasov_update(gkyl_vlasov_app* app, double dt)
{
  app->stat.nup += 1;
  struct timespec wst = gkyl_wall_clock();

  struct gkyl_update_status status = rk3(app, dt);
  app->tcurr += status.dt_actual;
  
  app->stat.total_tm += gkyl_time_diff_now_sec(wst);
  
  return status;
}

struct gkyl_vlasov_stat
gkyl_vlasov_app_stat(gkyl_vlasov_app* app)
{
  return app->stat;
}

void
gkyl_vlasov_app_species_ktm_rhs(gkyl_vlasov_app* app, int update_vol_term)
{
  for (int i=0; i<app->num_species; ++i) {
    
    struct vm_species *species = &app->species[i];
    
    const struct gkyl_array *qmem = app->field->qmem;
    gkyl_vlasov_set_qmem(species->eqn, qmem);

    const struct gkyl_array *fin = species->f;
    struct gkyl_array *rhs = species->f1;

    if (app->use_gpu)
      gkyl_hyper_dg_set_update_vol_cu(species->slvr, update_vol_term);
    else
      gkyl_hyper_dg_set_update_vol(species->slvr, update_vol_term);
    gkyl_array_clear_range(rhs, 0.0, species->local);
    if (app->use_gpu)
      gkyl_hyper_dg_advance_cu(species->slvr, species->local, fin,
        species->cflrate, rhs);
    else
      gkyl_hyper_dg_advance(species->slvr, species->local, fin,
        species->cflrate, rhs);
  }
}

void
gkyl_vlasov_app_stat_write(const gkyl_vlasov_app* app)
{
  const char *fmt = "%s-%s";
  int sz = gkyl_calc_strlen(fmt, app->name, "stat.json");
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, app->name, "stat.json");
  
  char buff[70];
  time_t t = time(NULL);
  struct tm curr_tm = *localtime(&t);

  // append to existing file so we have a history of different runs
  FILE *fp = 0;
  with_file (fp, fileNm, "a") {
    fprintf(fp, "{\n");

    if (strftime(buff, sizeof buff, "%c", &curr_tm))
      fprintf(fp, " \"date\" : \"%s\",\n", buff);

    fprintf(fp, " \"use_gpu\" : \"%d\",\n", app->stat.use_gpu);
    fprintf(fp, " \"nup\" : \"%ld\",\n", app->stat.nup);
    fprintf(fp, " \"nfeuler\" : \"%ld\",\n", app->stat.nfeuler);
    fprintf(fp, " \"nstage_2_fail\" : \"%ld\",\n", app->stat.nstage_2_fail);
    fprintf(fp, " \"nstage_3_fail\" : \"%ld\",\n", app->stat.nstage_3_fail);

    fprintf(fp, " \"stage_2_dt_diff\" : [ \"%lg\", \"%lg\" ],\n",
      app->stat.stage_2_dt_diff[0], app->stat.stage_2_dt_diff[1]);
    fprintf(fp, " \"stage_3_dt_diff\" : [ \"%lg\", \"%lg\" ],\n",
      app->stat.stage_3_dt_diff[0], app->stat.stage_3_dt_diff[1]);
    
    fprintf(fp, " \"total_tm\" : \"%lg\",\n", app->stat.total_tm);
    fprintf(fp, " \"init_species_tm\" : \"%lg\",\n", app->stat.init_species_tm);
    if (app->has_field)
      fprintf(fp, " \"init_field_tm\" : \"%lg\",\n", app->stat.init_field_tm);
    
    fprintf(fp, " \"species_rhs_tm\" : \"%lg\",\n", app->stat.species_rhs_tm);
    fprintf(fp, " \"species_coll_tm\" : \"%lg\",\n", app->stat.species_coll_tm);
    if (app->has_field) {
      fprintf(fp, " \"field_rhs_tm\" : \"%lg\",\n", app->stat.field_rhs_tm);
      fprintf(fp, " \"current_tm\" : \"%lg\",\n", app->stat.current_tm);
    }

    fprintf(fp, " \"nmom\" : \"%ld\",\n", app->stat.nmom);
    fprintf(fp, " \"mom_tm\" : \"%lg\"\n", app->stat.mom_tm);
  
    fprintf(fp, "}\n");
  }
}

void
gkyl_vlasov_app_release(gkyl_vlasov_app* app)
{
  gkyl_job_pool_release(app->job_pool);
  
  for (int i=0; i<app->num_species; ++i)
    vm_species_release(app, &app->species[i]);
  gkyl_free(app->species);
  if (app->has_field)
    vm_field_release(app, app->field);

  gkyl_free(app);
}
