#ifndef GKYL_ALLOC_H
#define GKYL_ALLOC_H

#include <stdlib.h>

#include <gkyl_util.h>

// The following allocators have the same calling/return behavior as
// standard C allocators. However, an error is signaled if allocation
// fails.

inline static void*
gkyl_malloc(size_t size)
{
  void *mem = malloc(size);
  if (NULL == mem) gkyl_exit("malloc failed!");
  return mem;
}

inline static void*
gkyl_calloc(size_t num, size_t size)
{
  void *mem = calloc(num, size);
  if (NULL == mem) gkyl_exit("calloc failed!");
  return mem;
}

inline static void*
gkyl_realloc(void *ptr, size_t new_size)
{
  void *mem = realloc(ptr, new_size);
  if (NULL == mem) gkyl_exit("realloc failed!");
  return mem;
}

inline static void
gkyl_free(void *ptr)
{
  free(ptr); ptr = 0;
}

/**
 * Allocate memory that is aligned to given byte boundary. You must
 * only use gkyl_aligned_realloc() and gkyl_aligned_free() methods for
 * memory returned by this method.
 *
 * @param align Alignment boundary. Must be power of 2.
 * @param size Number of bytes to allocate.
 */
void *gkyl_aligned_alloc(size_t align, size_t size);

/**
 * Reallocate memory that is aligned to given byte boundary. The
 * pointer passed to this must be allocated by
 * gkyl_aligned_alloc(). Returned pointer must be freed using
 * gkyl_aligned_free() method.
 *
 * @param align Alignment boundary. Must be power of 2.
 * @param size Number of bytes to allocate.
 */
//void *gkyl_aligned_realloc(size_t align, size_t size);

/**
 * Free memory allocated by gkyl_aligned_alloc().
 *
 * @param ptr Memory to free.
 */
void gkyl_aligned_free(void *ptr);

#endif // GKYL_ALLOC_H