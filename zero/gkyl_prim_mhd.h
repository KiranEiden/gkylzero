#pragma once

/**
 * Compute the fluid pressure given the conserved variables.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 */
static inline
double gkyl_mhd_pressure(double gas_gamma, const double q[8])
{
  return (gas_gamma-1) *
    (q[4] - 0.5*(q[1]*q[1]+q[2]*q[2]+q[3]*q[3])/q[0]
      - 0.5*(q[5]*q[5]+q[6]*q[6]+q[7]*q[7]));
}

/**
 * Compute maximum absolute speed.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param q Conserved variables
 * @return Maximum absolute speed for given q
 */
double gkyl_mhd_max_abs_speed(double gas_gamma, const double q[8]);

/**
 * Compute flux. Assumes rotation to local coordinate system.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
void gkyl_mhd_flux(double gas_gamma, const double q[8], double flux[8]);

/**
 * Compute flux for the GLM-MHD equations. Assumes rotation to local coordinate system.
 *
 * @param gas_gamma Gas adiabatic constant
 * @param ch The propagation speed of div(B) errors in the GLM method.
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
void gkyl_glm_mhd_flux(double gas_gamma, double ch, const double q[9], double flux[9]);
