#pragma once

/**
 * Compute maximum absolute speed.
 * 
 * @param c Speed of light
 * @param e_fact Correction speed for div(E) correction
 * @param b_fact Correction speed for div(B) correction
 * @param q Conserved variables
 * @return Maximum absolute speed for given q
 */
double gkyl_gr_maxwell_max_abs_speed(const double q[11]);

/**
 * Compute flux.
 * 
 * @param c Speed of light
 * @param e_fact Correction speed for div(E) correction
 * @param b_fact Correction speed for div(B) correction
 * @param Conserved variables
 * @param flux On output, the flux in direction 'dir'
 */
void gkyl_gr_maxwell_flux(const double q[11], double flux[11]);
