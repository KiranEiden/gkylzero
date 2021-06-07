#pragma once

#include <gkyl_ref_count.h>
#include <gkyl_util.h>

#include <stdbool.h>
#include <stdint.h>

// Type of element stored in array
enum gkyl_elem_type { GKYL_INT, GKYL_FLOAT, GKYL_DOUBLE, GKYL_USER };

/**
 * Array object. This is an untype, undimensioned, reference counted
 * array object. All additional structure is provided else where,
 * mainly by the range object.
 */
struct gkyl_array {
  enum gkyl_elem_type type; // type stored in array
  size_t elemsz, ncomp; // size of elements, number of 'components'
  size_t size; // number of indices
    
  uint32_t flags;
  size_t esznc; // elemsz*ncomp
  void *data; // pointer to data (do not use directly)
  struct gkyl_ref_count ref_count;

  // CUDA ONLY
  struct gkyl_array *on_device; // pointer to device clone if data is on device
  int numThreads, numBlocks; // CUDA kernel launch specifiers for entire array ops
};

/**
 * Create new array. Delete using gkyl_array_release method.
 * 
 * @param type Type of data in array
 * @param ncomp Number of components at each index
 * @param size Number of indices 
 * @return Pointer to newly allocated array.
 */
struct gkyl_array* gkyl_array_new(enum gkyl_elem_type type, size_t ncomp, size_t size);

/**
 * Create new array with data on NV-GPU. Delete using
 * gkyl_array_release method.
 *
 * NOTE: the data member lives on GPU, but the struct lives on the host.
 * For this reason, this method also creates a device struct and stores a 
 * pointer to it in the on_device member. This is a device clone of the 
 * host struct, and is what should be used to pass to CUDA kernels which
 * require the entire array struct on device.
 * 
 * @param type Type of data in array
 * @param ncomp Number of components at each index
 * @param size Number of indices 
 * @return Pointer to newly allocated array.
 */
struct gkyl_array* gkyl_array_cu_dev_new(enum gkyl_elem_type type, size_t ncomp, size_t size);

/**
 * Returns true if array lives on NV-GPU.
 *
 * @param arr Array to check
 * @return true of array lives on NV-GPU, false otherwise
 */
bool gkyl_array_is_cu_dev(const struct gkyl_array *arr);

/**
 * Copy into array: pointer to dest array is returned. 'dest' and
 * 'src' must not point to same data.
 * 
 * @param dest Destination for copy.
 * @param src Srouce to copy from.
 */
struct gkyl_array* gkyl_array_copy(struct gkyl_array* dest,
  const struct gkyl_array* src);

/**
 * Clone array: pointer to newly created array is returned.
 * 
 * @param arr Array to clone
 * @return Pointer to clone
 */
struct gkyl_array* gkyl_array_clone(const struct gkyl_array* arr);

/**
 * Fetches a pointer to the element stored at the index 'loc'.
 *
 * @param arr Array to fetch from
 * @param loc Element to fetch
 * @return Element at location 'loc'
 */
GKYL_CU_DH
static inline void*
gkyl_array_fetch(struct gkyl_array* arr, long loc)
{
  return ((char*) arr->data) + loc*arr->esznc;
}

/** Same as above, except fetches a constant pointer */
GKYL_CU_DH
static inline const void*
gkyl_array_cfetch(const struct gkyl_array* arr, long loc)
{
  return ((const char*) arr->data) + loc*arr->esznc;
}

/**
 * Aquire pointer to array. The pointer must be released using
 * gkyl_array_release method.
 *
 * @param arr Array to which a pointer is needed
 * @return Pointer to aquired array
 */
struct gkyl_array* gkyl_array_aquire(const struct gkyl_array* arr);

/**
 * Release pointer to array
 *
 * @param arr Array to release.
 */
void gkyl_array_release(const struct gkyl_array* arr);
