#pragma once

#include <gkyl_wv_eqn.h>

/**
 * Create a new Euler equation object.
 * 
 * @param gas_gamma Gas adiabatic constant
 * @return Pointer to Euler equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_euler_new(double gas_gamma);