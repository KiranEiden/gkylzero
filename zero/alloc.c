#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_util.h>

// Compute first 'align' boundary after 'num'
#define align_up(num, align)                    \
    (((num) + ((align) - 1)) & ~((align) - 1))

static const size_t PTR_OFFSET_SZ = sizeof(uint16_t);

void*
gkyl_malloc(size_t size)
{
  void *mem = malloc(size);
  if (0 == mem) gkyl_exit("malloc failed!");
  return mem;
}

void*
gkyl_calloc(size_t num, size_t size)
{
  void *mem = calloc(num, size);
  if (0 == mem) gkyl_exit("calloc failed!");
  return mem;
}

void*
gkyl_realloc(void *ptr, size_t new_size)
{
  void *mem = realloc(ptr, new_size);
  if (0 == mem) gkyl_exit("realloc failed!");
  return mem;
}

void
gkyl_free(void *ptr)
{
  free(ptr);
}

void*
gkyl_aligned_alloc(size_t align, size_t size)
{
  void *ptr = 0;
  assert((align & (align-1)) == 0); // power of 2?

  if (align && size) {
    uint32_t hdr_size = PTR_OFFSET_SZ + (align-1);
    void *p = gkyl_calloc(size+hdr_size, 1);
    if (p) {
      ptr = (void *) align_up(((uintptr_t)p + PTR_OFFSET_SZ), align);
      *((uint16_t *)ptr - 1) = (uint16_t) ((uintptr_t) ptr - (uintptr_t) p);
    }
  }
  return ptr;
}

void*
gkyl_aligned_realloc(void *ptr, size_t align, size_t old_sz, size_t new_sz)
{
  void *nptr = gkyl_aligned_alloc(align, new_sz);
  if (0 == nptr) {
    gkyl_exit("aligned_realloc failed!");
    return 0;
  }
  memcpy(nptr, ptr, old_sz < new_sz ? old_sz : new_sz);
  gkyl_aligned_free(ptr);
  return nptr;
}

void
gkyl_aligned_free(void* ptr)
{
  assert(ptr);
  uint16_t offset = *((uint16_t *)ptr - 1);
  gkyl_free((uint8_t *)ptr - offset);
}

// CUDA specific code

#ifdef GKYL_HAVE_CUDA

#include <cuda_runtime.h>

void*
gkyl_cu_malloc(size_t size)
{
  void *ptr;
  cudaError_t err = cudaMalloc(&ptr, size);
  if (err == cudaErrorMemoryAllocation)
    gkyl_exit("cudaMalloc failed!");
  return ptr;
}

void
gkyl_cu_free(void *ptr)
{
  cudaFree(ptr);
}

#endif // CUDA specific code
