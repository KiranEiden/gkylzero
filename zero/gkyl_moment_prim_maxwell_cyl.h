#pragma once

/**
 * Compute maximum absolute speed.
 * 
 * @param c Speed of light
 * @param q Conserved variables
 * @return Maximum absolute speed for given q
 */
double gkyl_maxwell_cyl_max_abs_speed(double c, const double q[9]);

/**
 * Compute flux.
 * 
 * @param c Speed of light
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
void gkyl_maxwell_cyl_flux(double c, const double q[9], double flux[9]);
