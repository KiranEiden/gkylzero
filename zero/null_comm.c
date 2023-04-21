#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_null_comm.h>

#include <string.h>

static void
comm_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_comm *comm = container_of(ref, struct gkyl_comm, ref_count);
  gkyl_free(comm);
}

static int
get_rank(struct gkyl_comm *comm, int *rank)
{
  *rank = 0;
  return 0;
}

static int
get_size(struct gkyl_comm *comm, int *sz)
{
  *sz = 1;
  return 0;
}

static int
all_reduce(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp,
  void *out)
{
  memcpy(out, inp, gkyl_elem_type_size[type]*nelem);
  return 0;
}

static int
array_sync(struct gkyl_comm *comm,
  const struct gkyl_range *local, const struct gkyl_range *local_ext,
  const int *nghost, struct gkyl_array *array)
{
  return 0;
}

static int barrier(struct gkyl_comm *comm) { return 0; }

static int
array_write(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  const struct gkyl_array *arr, const char *fname)
{
  return gkyl_grid_sub_array_write(grid, range, arr, fname);
}

struct gkyl_comm*
gkyl_null_comm_new(void)
{
  struct gkyl_comm *comm = gkyl_malloc(sizeof *comm);

  comm->get_rank = get_rank;
  comm->get_size = get_size;
  comm->all_reduce = all_reduce;
  comm->gkyl_array_sync = array_sync;
  comm->barrier = barrier;
  comm->gkyl_array_write = array_write;

  comm->ref_count = gkyl_ref_count_init(comm_free);

  return comm;
}