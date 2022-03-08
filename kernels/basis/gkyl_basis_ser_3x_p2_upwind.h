GKYL_CU_DH static inline void 
ser_3x_p2_upwind(const double* fUpwindQuad, double* GKYL_RESTRICT fUpwind) { 
  fUpwind[0] = 0.154320987654321*fUpwindQuad[8]+0.2469135802469136*fUpwindQuad[7]+0.154320987654321*fUpwindQuad[6]+0.2469135802469136*fUpwindQuad[5]+0.3950617283950617*fUpwindQuad[4]+0.2469135802469136*fUpwindQuad[3]+0.154320987654321*fUpwindQuad[2]+0.2469135802469136*fUpwindQuad[1]+0.154320987654321*fUpwindQuad[0]; 
  fUpwind[1] = 0.2070433312499806*fUpwindQuad[8]-0.2070433312499806*fUpwindQuad[6]+0.3312693299999688*fUpwindQuad[5]-0.3312693299999688*fUpwindQuad[3]+0.2070433312499806*fUpwindQuad[2]-0.2070433312499806*fUpwindQuad[0]; 
  fUpwind[2] = 0.2070433312499806*fUpwindQuad[8]+0.3312693299999688*fUpwindQuad[7]+0.2070433312499806*fUpwindQuad[6]-0.2070433312499806*fUpwindQuad[2]-0.3312693299999688*fUpwindQuad[1]-0.2070433312499806*fUpwindQuad[0]; 
  fUpwind[3] = 0.2777777777777778*fUpwindQuad[8]-0.2777777777777778*fUpwindQuad[6]-0.2777777777777778*fUpwindQuad[2]+0.2777777777777778*fUpwindQuad[0]; 
  fUpwind[4] = 0.138028887499987*fUpwindQuad[8]-0.2760577749999741*fUpwindQuad[7]+0.138028887499987*fUpwindQuad[6]+0.2208462199999792*fUpwindQuad[5]-0.4416924399999584*fUpwindQuad[4]+0.2208462199999792*fUpwindQuad[3]+0.138028887499987*fUpwindQuad[2]-0.2760577749999741*fUpwindQuad[1]+0.138028887499987*fUpwindQuad[0]; 
  fUpwind[5] = 0.138028887499987*fUpwindQuad[8]+0.2208462199999792*fUpwindQuad[7]+0.138028887499987*fUpwindQuad[6]-0.2760577749999741*fUpwindQuad[5]-0.4416924399999584*fUpwindQuad[4]-0.2760577749999741*fUpwindQuad[3]+0.138028887499987*fUpwindQuad[2]+0.2208462199999792*fUpwindQuad[1]+0.138028887499987*fUpwindQuad[0]; 
  fUpwind[6] = 0.1851851851851853*fUpwindQuad[8]-0.3703703703703705*fUpwindQuad[7]+0.1851851851851853*fUpwindQuad[6]-0.1851851851851853*fUpwindQuad[2]+0.3703703703703705*fUpwindQuad[1]-0.1851851851851853*fUpwindQuad[0]; 
  fUpwind[7] = 0.1851851851851853*fUpwindQuad[8]-0.1851851851851853*fUpwindQuad[6]-0.3703703703703705*fUpwindQuad[5]+0.3703703703703705*fUpwindQuad[3]+0.1851851851851853*fUpwindQuad[2]-0.1851851851851853*fUpwindQuad[0]; 
} 
