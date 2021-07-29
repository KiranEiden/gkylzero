#pragma once

#include <gkyl_util.h>

#include <stdio.h>
#include <stdlib.h>

// Represents a sized chunk of memory
struct gkyl_mem_chunk {
  size_t count; // size of memory in bytes
  char data[]; // Allocated memory
};

// The following allocators have the same calling/return behavior as
// standard C allocators. However, an error is signaled if allocation
// fails.

void* gkyl_malloc(size_t size);
void* gkyl_calloc(size_t num, size_t size);
void* gkyl_realloc(void *ptr, size_t new_size);
void gkyl_free(void *ptr);

/**
 * Allocate memory that is aligned to given byte boundary. You must
 * only use gkyl_aligned_realloc() and gkyl_aligned_free() methods to
 * reallocate or free memory returned by this method.
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
 * @param ptr Pointer to reallocate
 * @param align Alignment boundary. Must be power of 2.
 * @param old_sz Old size of memory.
 * @param new_sz New size of memory.
 */
void *gkyl_aligned_realloc(void *ptr, size_t align, size_t old_sz, size_t new_sz);

/**
 * Free memory allocated by gkyl_aligned_alloc().
 *
 * @param ptr Memory to free.
 */
void gkyl_aligned_free(void *ptr);

/**
 * Allocated a new chunk of count bytes of memory. Free using
 * gkyl_mem_chunk_release().
 *
 * @param count Bytes to allocte
 * @return Pointer to newly allocated chunk
 */
struct gkyl_mem_chunk* gkyl_mem_chunk_new(size_t count);

/**
 * Release memory chunk
 *
 * @param mem Memory chunk to release
 */
void gkyl_mem_chunk_release(struct gkyl_mem_chunk *mem);

// CUDA specific code (NV: Nvidia)

/** Allocate memory on NV-GPU */
void* gkyl_cu_malloc(size_t size);

/** Free memory on device */
void gkyl_cu_free(void *ptr);

/** Allocate pinned host memory on NV-GPU */
void* gkyl_cu_malloc_host(size_t size);

/** Free pinned host memory on device */
void gkyl_cu_free_host(void *ptr);

/** Copy data between host/device */
void gkyl_cu_memcpy(void *dst, void *src, size_t count, enum gkyl_cu_memcpy_kind kind);

/** Set memory on device */
void gkyl_cu_memset(void *data, int val, size_t count);
