#pragma once

// Identifiers for various equation systems
enum gkyl_eqn_type {
  GKYL_EQN_EULER, // Euler equations
  GKYL_EQN_SR_EULER, // SR Euler equations
  GKYL_EQN_ISO_EULER, // Isothermal Euler equations
  GKYL_EQN_TEN_MOMENT, // Ten-moment (with pressure tensor)
  GKYL_EQN_MAXWELL, // Maxwell equations
  GKYL_EQN_MHD,  // Ideal MHD equations
  GKYL_EQN_BURGERS, // Burgers equations
};

// Identifiers for specific field object types
enum gkyl_field_id {
  GKYL_FIELD_E_B = 0, // Maxwell (E, B). This is default
  GKYL_FIELD_SR_E_B, // Maxwell (E, B) with special relativity
  GKYL_FIELD_PHI, // Poisson (only phi)  
  GKYL_FIELD_PHI_A, // Poisson with static B = curl(A) (phi, A)
  GKYL_FIELD_NULL, // no field is present
  GKYL_FIELD_SR_NULL, // no field is present, special relativistic Vlasov
};

// Identifiers for specific diffusion object types
enum gkyl_diffusion_id {
  GKYL_ISO_DIFFUSION = 0, // Isotropic diffusion. This is default
  GKYL_ANISO_DIFFUSION, // Anisotropic diffusion. 
};

// Identifiers for specific collision object types
enum gkyl_collision_id {
  GKYL_NO_COLLISIONS = 0, // No collisions. This is default
  GKYL_BGK_COLLISIONS, // BGK Collision operator
  GKYL_LBO_COLLISIONS, // LBO Collision operator
  GKYL_FPO_COLLISIONS, // FPO Collision operator
};

// type of quadrature to use
enum gkyl_quad_type {
  GKYL_GAUSS_QUAD, // Gauss-Legendre quadrature
  GKYL_GAUSS_LOBATTO_QUAD, // Gauss-Lobatto quadrature
};

