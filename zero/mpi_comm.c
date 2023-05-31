
#ifdef GKYL_HAVE_MPI

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_format_desc.h>
#include <gkyl_comm_priv.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_mpi_comm.h>

#include <errno.h>

// Maximum number of recv neighbors: not sure hard-coding this is a
// good idea.
#define MAX_RECV_NEIGH 32

#define MPI_BASE_TAG 4242
#define MPI_BASE_PER_TAG 5252

// Mapping of Gkeyll type to MPI_Datatype
static MPI_Datatype g2_mpi_datatype[] = {
  [GKYL_INT] = MPI_INT,
  [GKYL_INT_64] = MPI_INT64_T,
  [GKYL_FLOAT] = MPI_FLOAT,
  [GKYL_DOUBLE] = MPI_DOUBLE
};

// Mapping of Gkeyll ops to MPI_Op
static MPI_Op g2_mpi_op[] = {
  [GKYL_MIN] = MPI_MIN,
  [GKYL_MAX] = MPI_MAX,
  [GKYL_SUM] = MPI_SUM
};

// Receive data
struct comm_buff_stat {
  struct gkyl_range range;
  MPI_Request status;
  gkyl_mem_buff buff;
};

// Private struct wrapping MPI-specific code
struct mpi_comm {
  struct gkyl_comm base; // base communicator

  MPI_Comm mcomm; // MPI communicator to use
  struct gkyl_rect_decomp *decomp; // pre-computed decomposition
  long local_range_offset; // offset of the local region

  struct gkyl_rect_decomp_neigh *neigh; // neighbors of local region
  struct gkyl_rect_decomp_neigh *per_neigh[GKYL_MAX_DIM]; // periodic neighbors

  struct gkyl_range dir_edge; // for use in computing tags
  int is_on_edge[2][GKYL_MAX_DIM]; // flags to indicate if local range is on edge
  bool touches_any_edge; // true if this range touches any edge

  int nrecv; // number of elements in rinfo array
  struct comm_buff_stat recv[MAX_RECV_NEIGH]; // info for recv data

  int nsend; // number of elements in sinfo array
  struct comm_buff_stat send[MAX_RECV_NEIGH]; // info for send data
};

static void
comm_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_comm *comm = container_of(ref, struct gkyl_comm, ref_count);
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, base);

  int ndim = mpi->decomp->ndim;
  gkyl_rect_decomp_release(mpi->decomp);

  gkyl_rect_decomp_neigh_release(mpi->neigh);
  for (int d=0; d<ndim; ++d)
    gkyl_rect_decomp_neigh_release(mpi->per_neigh[d]);

  for (int i=0; i<MAX_RECV_NEIGH; ++i)
    gkyl_mem_buff_release(mpi->recv[i].buff);

  for (int i=0; i<MAX_RECV_NEIGH; ++i)
    gkyl_mem_buff_release(mpi->send[i].buff);
  
  gkyl_free(mpi);
}

static int
get_rank(struct gkyl_comm *comm, int *rank)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, base);
  MPI_Comm_rank(mpi->mcomm, rank);
  return 0;
}

static int
get_size(struct gkyl_comm *comm, int *sz)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, base);
  MPI_Comm_size(mpi->mcomm, sz);
  return 0;
}

static int
all_reduce(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp,
  void *out)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, base);  
  int ret =
    MPI_Allreduce(inp, out, nelem, g2_mpi_datatype[type], g2_mpi_op[op], mpi->mcomm);
  return ret == MPI_SUCCESS ? 0 : 1;
}

static int
array_sync(struct gkyl_comm *comm,
  const struct gkyl_range *local, const struct gkyl_range *local_ext,
  struct gkyl_array *array)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, base);

  int rank;
  MPI_Comm_rank(mpi->mcomm, &rank);

  int elo[GKYL_MAX_DIM], eup[GKYL_MAX_DIM];
  for (int i=0; i<mpi->decomp->ndim; ++i)
    elo[i] = eup[i] = local_ext->upper[i]-local->upper[i];

  int nridx = 0;
  int tag = MPI_BASE_TAG;

    // post nonblocking recv to get data into ghost-cells  
  for (int n=0; n<mpi->neigh->num_neigh; ++n) {
    int nid = mpi->neigh->neigh[n];
    
    int isrecv = gkyl_sub_range_intersect(
      &mpi->recv[nridx].range, local_ext, &mpi->decomp->ranges[nid]);
    size_t recv_vol = array->esznc*mpi->recv[nridx].range.volume;

    if (isrecv) {
      if (gkyl_mem_buff_size(mpi->recv[nridx].buff) < recv_vol)
        gkyl_mem_buff_resize(mpi->recv[nridx].buff, recv_vol);
      
      MPI_Irecv(gkyl_mem_buff_data(mpi->recv[nridx].buff),
        recv_vol, MPI_CHAR, nid, tag, mpi->mcomm, &mpi->recv[nridx].status);

      nridx += 1;
    }
  }

  int nsidx = 0;
  
  // post non-blocking sends of skin-cell data to neighbors
  for (int n=0; n<mpi->neigh->num_neigh; ++n) {
    int nid = mpi->neigh->neigh[n];
    
    struct gkyl_range neigh_ext;
    gkyl_range_extend(&neigh_ext, &mpi->decomp->ranges[nid], elo, eup);

    int issend = gkyl_sub_range_intersect(
      &mpi->send[nsidx].range, local, &neigh_ext);
    size_t send_vol = array->esznc*mpi->send[nsidx].range.volume;

    if (issend) {
      if (gkyl_mem_buff_size(mpi->send[nsidx].buff) < send_vol)
        gkyl_mem_buff_resize(mpi->send[nsidx].buff, send_vol);
      
      gkyl_array_copy_to_buffer(gkyl_mem_buff_data(mpi->send[nsidx].buff),
        array, mpi->send[nsidx].range);
      
      MPI_Isend(gkyl_mem_buff_data(mpi->send[nsidx].buff),
        send_vol, MPI_CHAR, nid, tag, mpi->mcomm, &mpi->send[nsidx].status);

      nsidx += 1;
    }
  }

  // complete send
  for (int s=0; s<nsidx; ++s) {
    int issend = mpi->send[s].range.volume;
    if (issend)
      MPI_Wait(&mpi->send[s].status, MPI_STATUS_IGNORE);
  }

  // complete recv, copying data into ghost-cells
  for (int r=0; r<nridx; ++r) {
    int isrecv = mpi->recv[r].range.volume;
    if (isrecv) {
      MPI_Wait(&mpi->recv[r].status, MPI_STATUS_IGNORE);
      
      gkyl_array_copy_from_buffer(array,
        gkyl_mem_buff_data(mpi->recv[r].buff),
        mpi->recv[r].range
      );
    }
  }
  
  return 0;
}

static int
per_send_tag(const struct gkyl_range *dir_edge,
  int dir, int e)
{
  int base_tag = MPI_BASE_PER_TAG;
  return base_tag + gkyl_range_idx(dir_edge, (int[]) { dir, e });
}
static int
per_recv_tag(const struct gkyl_range *dir_edge,
  int dir, int e)
{
  int base_tag = MPI_BASE_PER_TAG;
  return base_tag + gkyl_range_idx(dir_edge, (int[]) { dir, (e+1)%2 });
}

static int
array_per_sync(struct gkyl_comm *comm, const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  int nper_dirs, const int *per_dirs, struct gkyl_array *array)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, base);

  if (!mpi->touches_any_edge) return 0; // nothing to sync

  int rank;
  MPI_Comm_rank(mpi->mcomm, &rank);
  
  int elo[GKYL_MAX_DIM], eup[GKYL_MAX_DIM];
  for (int i=0; i<mpi->decomp->ndim; ++i)
    elo[i] = eup[i] = local_ext->upper[i]-local->upper[i];

  int nridx = 0;
  int shift_sign[] = { -1, 1 };

  // post nonblocking recv to get data into ghost-cells
  for (int i=0; i<nper_dirs; ++i) {
    int dir = per_dirs[i];

    for (int e=0; e<2; ++e) {
      if (mpi->is_on_edge[e][dir]) {
        int delta[GKYL_MAX_DIM] = { 0 };
        delta[dir] = shift_sign[e]*gkyl_range_shape(&mpi->decomp->parent_range, dir);

        if (mpi->per_neigh[dir]->num_neigh == 1) { // really should  be a loop
          int nid = mpi->per_neigh[dir]->neigh[0];

          struct gkyl_range neigh_shift;
          gkyl_range_shift(&neigh_shift, &mpi->decomp->ranges[nid], delta);

          int isrecv = gkyl_sub_range_intersect(
            &mpi->recv[nridx].range, local_ext, &neigh_shift);
          size_t recv_vol = array->esznc*mpi->recv[nridx].range.volume;

          if (isrecv) {
            if (gkyl_mem_buff_size(mpi->recv[nridx].buff) < recv_vol)
              gkyl_mem_buff_resize(mpi->recv[nridx].buff, recv_vol);

            int tag = per_recv_tag(&mpi->dir_edge, dir, e);
            MPI_Irecv(gkyl_mem_buff_data(mpi->recv[nridx].buff),
              recv_vol, MPI_CHAR, nid, tag, mpi->mcomm, &mpi->recv[nridx].status);

            nridx += 1;
          }
        }
      }
    }
  }

  int nsidx = 0;  

  // post non-blocking send skin-cell data to neighbors
  for (int i=0; i<nper_dirs; ++i) {
    int dir = per_dirs[i];

    for (int e=0; e<2; ++e) {
      if (mpi->is_on_edge[e][dir]) {
        int delta[GKYL_MAX_DIM] = { 0 };
        delta[dir] = shift_sign[e]*gkyl_range_shape(&mpi->decomp->parent_range, dir);

        if (mpi->per_neigh[dir]->num_neigh == 1) { // really should  be a loop
          int nid = mpi->per_neigh[dir]->neigh[0];

          struct gkyl_range neigh_shift, neigh_shift_ext;
          gkyl_range_shift(&neigh_shift, &mpi->decomp->ranges[nid], delta);

          gkyl_range_extend(&neigh_shift_ext, &neigh_shift, elo, eup);
          int issend = gkyl_sub_range_intersect(
            &mpi->send[nsidx].range, local, &neigh_shift_ext);
          
          size_t send_vol = array->esznc*mpi->send[nsidx].range.volume;

          if (issend) {
            if (gkyl_mem_buff_size(mpi->send[nsidx].buff) < send_vol)
              gkyl_mem_buff_resize(mpi->send[nsidx].buff, send_vol);
            
            gkyl_array_copy_to_buffer(gkyl_mem_buff_data(mpi->send[nsidx].buff),
              array, mpi->send[nsidx].range);

            int tag = per_send_tag(&mpi->dir_edge, dir, e);
            MPI_Isend(gkyl_mem_buff_data(mpi->send[nsidx].buff),
              send_vol, MPI_CHAR, nid, tag, mpi->mcomm, &mpi->send[nsidx].status);

            nsidx += 1;
          }
        }
      }
    }
  }

  // complete send
  for (int s=0; s<nsidx; ++s) {
    int issend = mpi->send[s].range.volume;
    if (issend)
      MPI_Wait(&mpi->send[s].status, MPI_STATUS_IGNORE);
  }

  // complete recv, copying data into ghost-cells
  for (int r=0; r<nridx; ++r) {
    int isrecv = mpi->recv[r].range.volume;
    if (isrecv) {
      MPI_Wait(&mpi->recv[r].status, MPI_STATUS_IGNORE);
      
      gkyl_array_copy_from_buffer(array,
        gkyl_mem_buff_data(mpi->recv[r].buff),
        mpi->recv[r].range
      );
    }
  }

  mpi->nrecv = nridx > mpi->nrecv ? nridx : mpi->nrecv;
  
  return 0;
}

static int
barrier(struct gkyl_comm *comm)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, base);
  MPI_Barrier(mpi->mcomm);
  return 0;
}

// set of functions to help with parallel array output using MPI-IO
static void
sub_array_decomp_write(struct mpi_comm *comm, const struct gkyl_rect_decomp *decomp,
  const struct gkyl_range *range,
  const struct gkyl_array *arr, MPI_File fp)
{
#define _F(loc) gkyl_array_cfetch(arr, loc)

  int rank;
  MPI_Comm_rank(comm->mcomm, &rank);

  // seek to appropriate place in the file, depending on rank
  size_t hdr_sz = gkyl_base_hdr_size(0) + gkyl_file_type_3_hrd_size(range->ndim);
  size_t file_loc = hdr_sz +
    arr->esznc*comm->local_range_offset +
    rank*gkyl_file_type_3_range_hrd_size(range->ndim);

  MPI_Offset fp_offset = file_loc;
  MPI_File_seek(fp, fp_offset, MPI_SEEK_SET);

  do {
    char *buff; size_t buff_sz;
    FILE *fbuff = open_memstream(&buff, &buff_sz);
  
    uint64_t loidx[GKYL_MAX_DIM] = {0}, upidx[GKYL_MAX_DIM] = {0};
    for (int d = 0; d < range->ndim; ++d) {
      loidx[d] = range->lower[d];
      upidx[d] = range->upper[d];
    }
    
    fwrite(loidx, sizeof(uint64_t), range->ndim, fbuff);
    fwrite(upidx, sizeof(uint64_t), range->ndim, fbuff);
    uint64_t sz = range->volume;
    fwrite(&sz, sizeof(uint64_t), 1, fbuff);

    fclose(fbuff);

    MPI_Status status;
    MPI_File_write(fp, buff, buff_sz, MPI_CHAR, &status);

    free(buff);
    
  } while (0);

  // construct skip iterator to allow writing (potentially) in chunks
  // rather than element by element or requiring a copy of data. Note:
  // We must use "range" here and not decomp->ranges[rank] as the
  // latter is not a sub-range of the local extended
  // range.
  struct gkyl_range_skip_iter skip;
  gkyl_range_skip_iter_init(&skip, range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &skip.range);

  MPI_Status status;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&skip.range, iter.idx);
    MPI_File_write(fp, _F(start), arr->esznc*skip.delta, MPI_CHAR, &status);
  }

#undef _F
}

static int
grid_sub_array_decomp_write_fp(struct mpi_comm *comm,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_rect_decomp *decomp, const struct gkyl_range *range,
  const struct gkyl_array *arr, MPI_File fp)
{
  char *buff; size_t buff_sz;
  FILE *fbuff = open_memstream(&buff, &buff_sz);

  // wtite header to a char buffer
  gkyl_grid_sub_array_header_write_fp(grid,
    &(struct gkyl_array_header_info) {
      .file_type = gkyl_file_type_int[GKYL_MULTI_RANGE_DATA_FILE],
      .etype = arr->type,
      .esznc = arr->esznc,
      .tot_cells = decomp->parent_range.volume
    },
    fbuff
  );
  uint64_t nrange = decomp->ndecomp;
  fwrite(&nrange, sizeof(uint64_t), 1, fbuff);
  fclose(fbuff);

  // actual header IO is only done by rank 0
  int rank;
  MPI_Comm_rank(comm->mcomm, &rank);
  if (rank == 0) {
    MPI_Status status;
    MPI_File_write(fp, buff, buff_sz, MPI_CHAR, &status);
  }
  free(buff);
  
  // write data in array
  sub_array_decomp_write(comm, decomp, range, arr, fp);
  return errno;
}

static int
array_write(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  const struct gkyl_array *arr, const char *fname)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, base);
  MPI_File fp;
  int err =
    MPI_File_open(mpi->mcomm, fname, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
  if (err != MPI_SUCCESS)
    return err;
  err = grid_sub_array_decomp_write_fp(mpi, grid, mpi->decomp, range, arr, fp);
  MPI_File_close(&fp);
  return err;
}

static struct gkyl_comm*
extend_comm(const struct gkyl_comm *comm, const struct gkyl_range *erange)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, base);

  // extend internal decomp object and create a new communicator
  struct gkyl_rect_decomp *ext_decomp = gkyl_rect_decomp_extended_new(erange, mpi->decomp);
  struct gkyl_comm *ext_comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = mpi->mcomm,
      .decomp = ext_decomp
    }
  );
  gkyl_rect_decomp_release(ext_decomp);
  return ext_comm;
}

struct gkyl_comm*
gkyl_mpi_comm_new(const struct gkyl_mpi_comm_inp *inp)
{
  struct mpi_comm *mpi = gkyl_malloc(sizeof *mpi);
  mpi->mcomm = inp->mpi_comm;
  mpi->decomp = gkyl_rect_decomp_acquire(inp->decomp);

  int rank;
  MPI_Comm_rank(inp->mpi_comm, &rank);
  mpi->local_range_offset = gkyl_rect_decomp_calc_offset(mpi->decomp, rank);

  mpi->neigh = gkyl_rect_decomp_calc_neigh(mpi->decomp, inp->sync_corners, rank);
  for (int d=0; d<mpi->decomp->ndim; ++d)
    mpi->per_neigh[d] =
        gkyl_rect_decomp_calc_periodic_neigh(mpi->decomp, d, false, rank);

  gkyl_range_init(&mpi->dir_edge, 2, (int[]) { 0, 0 }, (int[]) { GKYL_MAX_DIM, 2 });

  int num_touches = 0;
  for (int d=0; d<mpi->decomp->ndim; ++d) {
    mpi->is_on_edge[0][d] = gkyl_range_is_on_lower_edge(
      d, &mpi->decomp->ranges[rank], &mpi->decomp->parent_range);
    mpi->is_on_edge[1][d] = gkyl_range_is_on_upper_edge(
      d, &mpi->decomp->ranges[rank], &mpi->decomp->parent_range);
    num_touches += mpi->is_on_edge[0][d] + mpi->is_on_edge[1][d];
  }
  mpi->touches_any_edge = num_touches > 0 ? true : false;
  
  mpi->nrecv = 0;
  for (int i=0; i<MAX_RECV_NEIGH; ++i)
    mpi->recv[i].buff = gkyl_mem_buff_new(16);

  mpi->nsend = 0;
  for (int i=0; i<MAX_RECV_NEIGH; ++i)
    mpi->send[i].buff = gkyl_mem_buff_new(16);
  
  mpi->base.get_rank = get_rank;
  mpi->base.get_size = get_size;
  mpi->base.barrier = barrier;
  mpi->base.all_reduce = all_reduce;
  mpi->base.gkyl_array_sync = array_sync;
  mpi->base.gkyl_array_per_sync = array_per_sync;
  mpi->base.gkyl_array_write = array_write;
  mpi->base.extend_comm = extend_comm;

  mpi->base.ref_count = gkyl_ref_count_init(comm_free);

  return &mpi->base;
}

#endif
