#include <float.h>

#include <acutest.h>
#include <gkyl_update_fsm.h>

enum seq_states {
  SRC_1 = GKYL_UPDATE_FSM_FIRST, // first step
  FLUID,
  SRC_2,
  FINIS // sentinel
};

struct seq_ctx {
  int nredo, nfluid, nsrc;
};

struct gkyl_update_status
seq_redo(double tcurr, double dt, void *ctx)
{
  //printf("seq_redo: dt = %g\n", dt);
  
  struct seq_ctx *sc = ctx;
  sc->nredo += 1;
  
  return (struct gkyl_update_status) {
    .next_state = SRC_1,
    .dt_actual = dt,
    .dt_suggested = DBL_MAX
  };
}

struct gkyl_update_status
seq_src_1(double tcurr, double dt, void *ctx)
{
  //printf("seq_src_1: dt = %g\n", dt);
  
  struct seq_ctx *sc = ctx;
  sc->nsrc += 1;

  int next_state = FLUID;
  if (dt > 100)
    next_state = GKYL_UPDATE_FSM_ABORT;  

  return (struct gkyl_update_status) {
    .next_state = next_state,
    .dt_actual = dt,
    .dt_suggested = DBL_MAX,
  };
}

struct gkyl_update_status
seq_fluid(double tcurr, double dt, void *ctx)
{
  //printf("seq_fluid: dt = %g\n", dt);
  
  struct seq_ctx *sc = ctx;
  sc->nfluid += 1;

  double max_dt = 0.1, dt_actual = dt;

  if (dt > max_dt)
    dt_actual = max_dt;

  // take time-step of dt_actual
  
  return (struct gkyl_update_status) {
    .next_state = SRC_2,
    .dt_actual = dt_actual,
    .dt_suggested = max_dt,
  };
}

struct gkyl_update_status
seq_src_2(double tcurr, double dt, void *ctx)
{
  //printf("seq_src_2: dt = %g\n", dt);
  
  struct seq_ctx *sc = ctx;
  sc->nsrc += 1;
  
  return (struct gkyl_update_status) {
    .next_state = GKYL_UPDATE_FSM_FINISH,
    .dt_actual = dt,
    .dt_suggested = DBL_MAX,
  };
}

void
test_seq_1()
{
  struct seq_ctx ctx = { };

  struct gkyl_update_fsm *seq = gkyl_update_fsm_new(
    FINIS, // number of steps in sequence
    (struct gkyl_update_fsm_step[FINIS]) { // steps in sequence
      [GKYL_UPDATE_FSM_REDO] = { .ctx = &ctx, .u = seq_redo },
      [SRC_1] = { .ctx = &ctx, .u = seq_src_1 },
      [FLUID] = { .ctx = &ctx, .u = seq_fluid },
      [SRC_2] = { .ctx = &ctx, .u = seq_src_2 }
    }
  );

  TEST_CHECK( 4 == seq->nsteps );
  
  struct gkyl_update_status status = gkyl_update_fsm_run(seq, 0.0, 1.0);

  TEST_CHECK( 1 == ctx.nredo );
  TEST_CHECK( 2 == ctx.nfluid );
  TEST_CHECK( 3 == ctx.nsrc );

  TEST_CHECK( 0.1 == status.dt_actual );
  TEST_CHECK( 0.1 == status.dt_suggested );

  status = gkyl_update_fsm_run(seq, 0.0, 200.0); // should abort
  TEST_CHECK( false == status.success );

  gkyl_update_fsm_release(seq);
}

TEST_LIST = {
  { "seq_1", test_seq_1 },
  { NULL, NULL },  
};
