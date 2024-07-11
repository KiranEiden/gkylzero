#pragma once

#include <gkyl_array.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// Spatial coordinate system flag
enum gkyl_wave_coord_flag {
  WAVE_COORD_CART = 0, // Cartesian (default)
  WAVE_COORD_CYL, // Cylindrical (r, phi, z) with unit basis; 2D polar if ndim = 2
  WAVE_COORD_SPH // Spherical (r, theta, phi) with unit basis
};

// Functions for handling different coordinate systems and spacetimes
struct gkyl_wave_geom_inp {
  evalf_t mapc2p; // Map computational to physical coordinates
  evalf3_t get_cov_basis; // Calculate covariant computational basis vectors
  evalf3_t get_con_basis; // Calculate contravariant computational basis vectors
  struct gkyl_gr_spacetime *spacetime; // Pointer to spacetime struct (if provided, others will be ignored)
};

// Geometry information for a single cell: recall a cell "owns" the
// faces on the lower side of cell. We store both covariant and contravariant vector components
struct gkyl_wave_cell_geom {
  double kappa; // ratio of cell-volume in phy to comp space
  double lenr[GKYL_MAX_CDIM]; // ratio of face-area in phys to comp space for "lower" faces
  double chr_sym[GKYL_MAX_CDIM][GKYL_MAX_CDIM][GKYL_MAX_CDIM]; // Spatial Christoffel symbols with last index first 
  // norm[d] is the normal to face perp to direction 'd'
  double norm_cov[GKYL_MAX_CDIM][GKYL_MAX_CDIM]; // Covariant components
  double norm_con[GKYL_MAX_CDIM][GKYL_MAX_CDIM]; // Contravariant components
  // tau1[d] X tau2[d] = norm[d] are tangents to face perp to direction 'd'
  double tau1_cov[GKYL_MAX_CDIM][GKYL_MAX_CDIM];
  double tau1_con[GKYL_MAX_CDIM][GKYL_MAX_CDIM];
  double tau2_cov[GKYL_MAX_CDIM][GKYL_MAX_CDIM];
  double tau2_con[GKYL_MAX_CDIM][GKYL_MAX_CDIM];
};

// geometry information over a range of cells
struct gkyl_wave_geom {
  struct gkyl_range range; // range over which geometry is defined
  struct gkyl_array *geom; // geometry in each cell

  uint32_t flags;
  struct gkyl_ref_count ref_count;  
  struct gkyl_wave_geom *on_dev; // pointer to itself or device object
};

/**
 * Create a new wave geometry object. 
 *
 * @param grid Grid on which geometry lives
 * @param range Range on which geometry should be constructed
 * @param mapc2p Mapping from computational to physical space
 * @param ctx Context for use in mapping
 */
struct gkyl_wave_geom*
gkyl_wave_geom_new(const struct gkyl_rect_grid *grid,
  struct gkyl_range *range, evalf_t mapc2p, void *ctx, bool use_gpu);

/**
 * Set components of coordinate mapping struct based on predefined coordinate flag.
 *
 * @param cflag Flag for predefined coordinate system
 * @param inp Pointer to input struct to set up
 */
void
gkyl_wave_geom_inp_from_flag(enum gkyl_wave_coord_flag cflag, struct gkyl_wave_geom_inp *inp);

/**
 * Create a new wave geometry object from coordinate system flag.
 *
 * @param grid Grid on which geometry lives
 * @param range Range on which geometry should be constructed
 * @param cflag Flag for predefined coordinate system -- sets mapc2p and basis vectors
 * @param ctx Context for use in mapping
 */
struct gkyl_wave_geom*
gkyl_wave_geom_from_coord_flag(const struct gkyl_rect_grid *grid,
  struct gkyl_range *range, enum gkyl_wave_coord_flag cflag, void *ctx, bool use_gpu);
  
/**
 * Create a new wave geometry object from coordinate map struct.
 *
 * @param grid Grid on which geometry lives
 * @param range Range on which geometry should be constructed
 * @param inp Struct containing either mapc2p and basis vector calculation functions or a spacetime struct.
 * @param ctx Context for use in mapping
 */
struct gkyl_wave_geom*
gkyl_wave_geom_from_wg_inp(const struct gkyl_rect_grid *grid,
  struct gkyl_range *range, struct gkyl_wave_geom_inp *inp, void *ctx, bool use_gpu);

/**
 * Create a new wave geometry object that lives on NV-GPU: see new() method
 * above for documentation.
 */
struct gkyl_wave_geom*
gkyl_wave_geom_cu_dev_new(const struct gkyl_rect_grid *grid,
  struct gkyl_range *range, evalf_t mapc2p, void *ctx);

/**
 * Create a new wave geometry object that lives on NV-GPU: see ...from_coord_flag() method
 * above for documentation.
 */
 struct gkyl_wave_geom*
gkyl_wave_geom_cu_dev_from_coord_flag(const struct gkyl_rect_grid *grid,
  struct gkyl_range *range, enum gkyl_wave_coord_flag cflag, void *ctx);
  
/**
 * Create a new wave geometry object that lives on NV-GPU: see ...from_coord_maps() method
 * above for documentation.
 */
struct gkyl_wave_geom*
gkyl_wave_geom_cu_dev_from_wg_inp(const struct gkyl_rect_grid *grid,
  struct gkyl_range *range, struct gkyl_wave_geom_inp *inp, void *ctx);

/**
 * Acquire pointer to geometry object. The pointer must be released
 * using gkyl_wave_geom_release method.
 *
 * @param wg Geometry to which a pointer is needed
 * @return Pointer to acquired geometry
 */
struct gkyl_wave_geom* gkyl_wave_geom_acquire(const struct gkyl_wave_geom* wg);

/**
 * Get pointer to geometry in a cell given by idx into the range over
 * which the geometry was constructed.
 */
GKYL_CU_DH
static inline const struct gkyl_wave_cell_geom*
gkyl_wave_geom_get(const struct gkyl_wave_geom *wg, const int *idx)
{
  return (const struct gkyl_wave_cell_geom*) gkyl_array_cfetch(wg->geom, gkyl_range_idx(&wg->range, idx));
}

/**
 * Release geometry object.
 *
 * @param wg Wave geometry object to release.
 */
void gkyl_wave_geom_release(const struct gkyl_wave_geom *wg);
