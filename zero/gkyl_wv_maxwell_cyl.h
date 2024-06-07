#pragma once

#include <gkyl_wv_eqn.h>

/**
 * Create a new maxwell equation object.
 * 
 * @param c_speed Speed of light
 * @param e_fact Factor of light-speed for electric field correction
 * @param b_fact Factor of light-speed for magnetic field correction
 * @return Pointer to Maxwell equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_maxwell_cyl_new(double c);
