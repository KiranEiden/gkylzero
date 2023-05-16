#include <gkyl_fem_poisson_perp_kernels.h> 
 
GKYL_CU_DH void fem_poisson_perp_sol_stencil_3x_ser_p1(const double *sol_nodal_global, const long *globalIdxs, double *sol_modal_local) 
{ 
  // sol_nodal_global: global nodal solution vector.
  // sol_modal_local: local modal solution vector.

  sol_modal_local[0] = 0.3535533905932737*sol_nodal_global[globalIdxs[7]]+0.3535533905932737*sol_nodal_global[globalIdxs[6]]+0.3535533905932737*sol_nodal_global[globalIdxs[5]]+0.3535533905932737*sol_nodal_global[globalIdxs[4]]+0.3535533905932737*sol_nodal_global[globalIdxs[3]]+0.3535533905932737*sol_nodal_global[globalIdxs[2]]+0.3535533905932737*sol_nodal_global[globalIdxs[1]]+0.3535533905932737*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[1] = 0.2041241452319315*sol_nodal_global[globalIdxs[7]]-0.2041241452319315*sol_nodal_global[globalIdxs[6]]+0.2041241452319315*sol_nodal_global[globalIdxs[5]]-0.2041241452319315*sol_nodal_global[globalIdxs[4]]+0.2041241452319315*sol_nodal_global[globalIdxs[3]]-0.2041241452319315*sol_nodal_global[globalIdxs[2]]+0.2041241452319315*sol_nodal_global[globalIdxs[1]]-0.2041241452319315*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[2] = 0.2041241452319315*sol_nodal_global[globalIdxs[7]]+0.2041241452319315*sol_nodal_global[globalIdxs[6]]-0.2041241452319315*sol_nodal_global[globalIdxs[5]]-0.2041241452319315*sol_nodal_global[globalIdxs[4]]+0.2041241452319315*sol_nodal_global[globalIdxs[3]]+0.2041241452319315*sol_nodal_global[globalIdxs[2]]-0.2041241452319315*sol_nodal_global[globalIdxs[1]]-0.2041241452319315*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[3] = 0.2041241452319315*sol_nodal_global[globalIdxs[7]]+0.2041241452319315*sol_nodal_global[globalIdxs[6]]+0.2041241452319315*sol_nodal_global[globalIdxs[5]]+0.2041241452319315*sol_nodal_global[globalIdxs[4]]-0.2041241452319315*sol_nodal_global[globalIdxs[3]]-0.2041241452319315*sol_nodal_global[globalIdxs[2]]-0.2041241452319315*sol_nodal_global[globalIdxs[1]]-0.2041241452319315*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[4] = 0.1178511301977579*sol_nodal_global[globalIdxs[7]]-0.1178511301977579*sol_nodal_global[globalIdxs[6]]-0.1178511301977579*sol_nodal_global[globalIdxs[5]]+0.1178511301977579*sol_nodal_global[globalIdxs[4]]+0.1178511301977579*sol_nodal_global[globalIdxs[3]]-0.1178511301977579*sol_nodal_global[globalIdxs[2]]-0.1178511301977579*sol_nodal_global[globalIdxs[1]]+0.1178511301977579*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[5] = 0.1178511301977579*sol_nodal_global[globalIdxs[7]]-0.1178511301977579*sol_nodal_global[globalIdxs[6]]+0.1178511301977579*sol_nodal_global[globalIdxs[5]]-0.1178511301977579*sol_nodal_global[globalIdxs[4]]-0.1178511301977579*sol_nodal_global[globalIdxs[3]]+0.1178511301977579*sol_nodal_global[globalIdxs[2]]-0.1178511301977579*sol_nodal_global[globalIdxs[1]]+0.1178511301977579*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[6] = 0.1178511301977579*sol_nodal_global[globalIdxs[7]]+0.1178511301977579*sol_nodal_global[globalIdxs[6]]-0.1178511301977579*sol_nodal_global[globalIdxs[5]]-0.1178511301977579*sol_nodal_global[globalIdxs[4]]-0.1178511301977579*sol_nodal_global[globalIdxs[3]]-0.1178511301977579*sol_nodal_global[globalIdxs[2]]+0.1178511301977579*sol_nodal_global[globalIdxs[1]]+0.1178511301977579*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[7] = 0.06804138174397716*sol_nodal_global[globalIdxs[7]]-0.06804138174397716*sol_nodal_global[globalIdxs[6]]-0.06804138174397716*sol_nodal_global[globalIdxs[5]]+0.06804138174397716*sol_nodal_global[globalIdxs[4]]-0.06804138174397716*sol_nodal_global[globalIdxs[3]]+0.06804138174397716*sol_nodal_global[globalIdxs[2]]+0.06804138174397716*sol_nodal_global[globalIdxs[1]]-0.06804138174397716*sol_nodal_global[globalIdxs[0]];

}
GKYL_CU_DH void fem_poisson_perp_sol_stencil_3x_ser_p2(const double *sol_nodal_global, const long *globalIdxs, double *sol_modal_local) 
{ 
  // sol_nodal_global: global nodal solution vector.
  // sol_modal_local: local modal solution vector.

  sol_modal_local[0] = (-0.3535533905932737*sol_nodal_global[globalIdxs[19]])+0.4714045207910316*sol_nodal_global[globalIdxs[18]]-0.3535533905932737*sol_nodal_global[globalIdxs[17]]+0.4714045207910316*sol_nodal_global[globalIdxs[16]]+0.4714045207910316*sol_nodal_global[globalIdxs[15]]-0.3535533905932737*sol_nodal_global[globalIdxs[14]]+0.4714045207910316*sol_nodal_global[globalIdxs[13]]-0.3535533905932737*sol_nodal_global[globalIdxs[12]]+0.4714045207910316*sol_nodal_global[globalIdxs[11]]+0.4714045207910316*sol_nodal_global[globalIdxs[10]]+0.4714045207910316*sol_nodal_global[globalIdxs[9]]+0.4714045207910316*sol_nodal_global[globalIdxs[8]]-0.3535533905932737*sol_nodal_global[globalIdxs[7]]+0.4714045207910316*sol_nodal_global[globalIdxs[6]]-0.3535533905932737*sol_nodal_global[globalIdxs[5]]+0.4714045207910316*sol_nodal_global[globalIdxs[4]]+0.4714045207910316*sol_nodal_global[globalIdxs[3]]-0.3535533905932737*sol_nodal_global[globalIdxs[2]]+0.4714045207910316*sol_nodal_global[globalIdxs[1]]-0.3535533905932737*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[1] = (-0.06804138174397717*sol_nodal_global[globalIdxs[19]])+0.06804138174397717*sol_nodal_global[globalIdxs[17]]+0.2721655269759087*sol_nodal_global[globalIdxs[16]]-0.2721655269759087*sol_nodal_global[globalIdxs[15]]-0.06804138174397717*sol_nodal_global[globalIdxs[14]]+0.06804138174397717*sol_nodal_global[globalIdxs[12]]+0.2721655269759087*sol_nodal_global[globalIdxs[11]]-0.2721655269759087*sol_nodal_global[globalIdxs[10]]+0.2721655269759087*sol_nodal_global[globalIdxs[9]]-0.2721655269759087*sol_nodal_global[globalIdxs[8]]-0.06804138174397717*sol_nodal_global[globalIdxs[7]]+0.06804138174397717*sol_nodal_global[globalIdxs[5]]+0.2721655269759087*sol_nodal_global[globalIdxs[4]]-0.2721655269759087*sol_nodal_global[globalIdxs[3]]-0.06804138174397717*sol_nodal_global[globalIdxs[2]]+0.06804138174397717*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[2] = (-0.06804138174397717*sol_nodal_global[globalIdxs[19]])+0.2721655269759087*sol_nodal_global[globalIdxs[18]]-0.06804138174397717*sol_nodal_global[globalIdxs[17]]+0.06804138174397717*sol_nodal_global[globalIdxs[14]]-0.2721655269759087*sol_nodal_global[globalIdxs[13]]+0.06804138174397717*sol_nodal_global[globalIdxs[12]]+0.2721655269759087*sol_nodal_global[globalIdxs[11]]+0.2721655269759087*sol_nodal_global[globalIdxs[10]]-0.2721655269759087*sol_nodal_global[globalIdxs[9]]-0.2721655269759087*sol_nodal_global[globalIdxs[8]]-0.06804138174397717*sol_nodal_global[globalIdxs[7]]+0.2721655269759087*sol_nodal_global[globalIdxs[6]]-0.06804138174397717*sol_nodal_global[globalIdxs[5]]+0.06804138174397717*sol_nodal_global[globalIdxs[2]]-0.2721655269759087*sol_nodal_global[globalIdxs[1]]+0.06804138174397717*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[3] = (-0.06804138174397717*sol_nodal_global[globalIdxs[19]])+0.2721655269759087*sol_nodal_global[globalIdxs[18]]-0.06804138174397717*sol_nodal_global[globalIdxs[17]]+0.2721655269759087*sol_nodal_global[globalIdxs[16]]+0.2721655269759087*sol_nodal_global[globalIdxs[15]]-0.06804138174397717*sol_nodal_global[globalIdxs[14]]+0.2721655269759087*sol_nodal_global[globalIdxs[13]]-0.06804138174397717*sol_nodal_global[globalIdxs[12]]+0.06804138174397717*sol_nodal_global[globalIdxs[7]]-0.2721655269759087*sol_nodal_global[globalIdxs[6]]+0.06804138174397717*sol_nodal_global[globalIdxs[5]]-0.2721655269759087*sol_nodal_global[globalIdxs[4]]-0.2721655269759087*sol_nodal_global[globalIdxs[3]]+0.06804138174397717*sol_nodal_global[globalIdxs[2]]-0.2721655269759087*sol_nodal_global[globalIdxs[1]]+0.06804138174397717*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[4] = 0.0392837100659193*sol_nodal_global[globalIdxs[19]]-0.0392837100659193*sol_nodal_global[globalIdxs[17]]-0.0392837100659193*sol_nodal_global[globalIdxs[14]]+0.0392837100659193*sol_nodal_global[globalIdxs[12]]+0.1571348402636772*sol_nodal_global[globalIdxs[11]]-0.1571348402636772*sol_nodal_global[globalIdxs[10]]-0.1571348402636772*sol_nodal_global[globalIdxs[9]]+0.1571348402636772*sol_nodal_global[globalIdxs[8]]+0.0392837100659193*sol_nodal_global[globalIdxs[7]]-0.0392837100659193*sol_nodal_global[globalIdxs[5]]-0.0392837100659193*sol_nodal_global[globalIdxs[2]]+0.0392837100659193*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[5] = 0.0392837100659193*sol_nodal_global[globalIdxs[19]]-0.0392837100659193*sol_nodal_global[globalIdxs[17]]+0.1571348402636772*sol_nodal_global[globalIdxs[16]]-0.1571348402636772*sol_nodal_global[globalIdxs[15]]+0.0392837100659193*sol_nodal_global[globalIdxs[14]]-0.0392837100659193*sol_nodal_global[globalIdxs[12]]-0.0392837100659193*sol_nodal_global[globalIdxs[7]]+0.0392837100659193*sol_nodal_global[globalIdxs[5]]-0.1571348402636772*sol_nodal_global[globalIdxs[4]]+0.1571348402636772*sol_nodal_global[globalIdxs[3]]-0.0392837100659193*sol_nodal_global[globalIdxs[2]]+0.0392837100659193*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[6] = 0.0392837100659193*sol_nodal_global[globalIdxs[19]]+0.1571348402636772*sol_nodal_global[globalIdxs[18]]+0.0392837100659193*sol_nodal_global[globalIdxs[17]]-0.0392837100659193*sol_nodal_global[globalIdxs[14]]-0.1571348402636772*sol_nodal_global[globalIdxs[13]]-0.0392837100659193*sol_nodal_global[globalIdxs[12]]-0.0392837100659193*sol_nodal_global[globalIdxs[7]]-0.1571348402636772*sol_nodal_global[globalIdxs[6]]-0.0392837100659193*sol_nodal_global[globalIdxs[5]]+0.0392837100659193*sol_nodal_global[globalIdxs[2]]+0.1571348402636772*sol_nodal_global[globalIdxs[1]]+0.0392837100659193*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[7] = 0.105409255338946*sol_nodal_global[globalIdxs[19]]-0.2108185106778919*sol_nodal_global[globalIdxs[18]]+0.105409255338946*sol_nodal_global[globalIdxs[17]]+0.105409255338946*sol_nodal_global[globalIdxs[14]]-0.2108185106778919*sol_nodal_global[globalIdxs[13]]+0.105409255338946*sol_nodal_global[globalIdxs[12]]+0.105409255338946*sol_nodal_global[globalIdxs[7]]-0.2108185106778919*sol_nodal_global[globalIdxs[6]]+0.105409255338946*sol_nodal_global[globalIdxs[5]]+0.105409255338946*sol_nodal_global[globalIdxs[2]]-0.2108185106778919*sol_nodal_global[globalIdxs[1]]+0.105409255338946*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[8] = 0.105409255338946*sol_nodal_global[globalIdxs[19]]+0.105409255338946*sol_nodal_global[globalIdxs[17]]-0.2108185106778919*sol_nodal_global[globalIdxs[16]]-0.2108185106778919*sol_nodal_global[globalIdxs[15]]+0.105409255338946*sol_nodal_global[globalIdxs[14]]+0.105409255338946*sol_nodal_global[globalIdxs[12]]+0.105409255338946*sol_nodal_global[globalIdxs[7]]+0.105409255338946*sol_nodal_global[globalIdxs[5]]-0.2108185106778919*sol_nodal_global[globalIdxs[4]]-0.2108185106778919*sol_nodal_global[globalIdxs[3]]+0.105409255338946*sol_nodal_global[globalIdxs[2]]+0.105409255338946*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[9] = 0.105409255338946*sol_nodal_global[globalIdxs[19]]+0.105409255338946*sol_nodal_global[globalIdxs[17]]+0.105409255338946*sol_nodal_global[globalIdxs[14]]+0.105409255338946*sol_nodal_global[globalIdxs[12]]-0.2108185106778919*sol_nodal_global[globalIdxs[11]]-0.2108185106778919*sol_nodal_global[globalIdxs[10]]-0.2108185106778919*sol_nodal_global[globalIdxs[9]]-0.2108185106778919*sol_nodal_global[globalIdxs[8]]+0.105409255338946*sol_nodal_global[globalIdxs[7]]+0.105409255338946*sol_nodal_global[globalIdxs[5]]+0.105409255338946*sol_nodal_global[globalIdxs[2]]+0.105409255338946*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[10] = 0.06804138174397717*sol_nodal_global[globalIdxs[19]]-0.06804138174397717*sol_nodal_global[globalIdxs[17]]-0.06804138174397717*sol_nodal_global[globalIdxs[14]]+0.06804138174397717*sol_nodal_global[globalIdxs[12]]-0.06804138174397717*sol_nodal_global[globalIdxs[7]]+0.06804138174397717*sol_nodal_global[globalIdxs[5]]+0.06804138174397717*sol_nodal_global[globalIdxs[2]]-0.06804138174397717*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[11] = 0.06085806194501844*sol_nodal_global[globalIdxs[19]]-0.1217161238900369*sol_nodal_global[globalIdxs[18]]+0.06085806194501844*sol_nodal_global[globalIdxs[17]]-0.06085806194501844*sol_nodal_global[globalIdxs[14]]+0.1217161238900369*sol_nodal_global[globalIdxs[13]]-0.06085806194501844*sol_nodal_global[globalIdxs[12]]+0.06085806194501844*sol_nodal_global[globalIdxs[7]]-0.1217161238900369*sol_nodal_global[globalIdxs[6]]+0.06085806194501844*sol_nodal_global[globalIdxs[5]]-0.06085806194501844*sol_nodal_global[globalIdxs[2]]+0.1217161238900369*sol_nodal_global[globalIdxs[1]]-0.06085806194501844*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[12] = 0.06085806194501844*sol_nodal_global[globalIdxs[19]]-0.06085806194501844*sol_nodal_global[globalIdxs[17]]-0.1217161238900369*sol_nodal_global[globalIdxs[16]]+0.1217161238900369*sol_nodal_global[globalIdxs[15]]+0.06085806194501844*sol_nodal_global[globalIdxs[14]]-0.06085806194501844*sol_nodal_global[globalIdxs[12]]+0.06085806194501844*sol_nodal_global[globalIdxs[7]]-0.06085806194501844*sol_nodal_global[globalIdxs[5]]-0.1217161238900369*sol_nodal_global[globalIdxs[4]]+0.1217161238900369*sol_nodal_global[globalIdxs[3]]+0.06085806194501844*sol_nodal_global[globalIdxs[2]]-0.06085806194501844*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[13] = 0.06085806194501844*sol_nodal_global[globalIdxs[19]]-0.1217161238900369*sol_nodal_global[globalIdxs[18]]+0.06085806194501844*sol_nodal_global[globalIdxs[17]]+0.06085806194501844*sol_nodal_global[globalIdxs[14]]-0.1217161238900369*sol_nodal_global[globalIdxs[13]]+0.06085806194501844*sol_nodal_global[globalIdxs[12]]-0.06085806194501844*sol_nodal_global[globalIdxs[7]]+0.1217161238900369*sol_nodal_global[globalIdxs[6]]-0.06085806194501844*sol_nodal_global[globalIdxs[5]]-0.06085806194501844*sol_nodal_global[globalIdxs[2]]+0.1217161238900369*sol_nodal_global[globalIdxs[1]]-0.06085806194501844*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[14] = 0.06085806194501844*sol_nodal_global[globalIdxs[19]]+0.06085806194501844*sol_nodal_global[globalIdxs[17]]-0.1217161238900369*sol_nodal_global[globalIdxs[16]]-0.1217161238900369*sol_nodal_global[globalIdxs[15]]+0.06085806194501844*sol_nodal_global[globalIdxs[14]]+0.06085806194501844*sol_nodal_global[globalIdxs[12]]-0.06085806194501844*sol_nodal_global[globalIdxs[7]]-0.06085806194501844*sol_nodal_global[globalIdxs[5]]+0.1217161238900369*sol_nodal_global[globalIdxs[4]]+0.1217161238900369*sol_nodal_global[globalIdxs[3]]-0.06085806194501844*sol_nodal_global[globalIdxs[2]]-0.06085806194501844*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[15] = 0.06085806194501844*sol_nodal_global[globalIdxs[19]]-0.06085806194501844*sol_nodal_global[globalIdxs[17]]+0.06085806194501844*sol_nodal_global[globalIdxs[14]]-0.06085806194501844*sol_nodal_global[globalIdxs[12]]-0.1217161238900369*sol_nodal_global[globalIdxs[11]]+0.1217161238900369*sol_nodal_global[globalIdxs[10]]-0.1217161238900369*sol_nodal_global[globalIdxs[9]]+0.1217161238900369*sol_nodal_global[globalIdxs[8]]+0.06085806194501844*sol_nodal_global[globalIdxs[7]]-0.06085806194501844*sol_nodal_global[globalIdxs[5]]+0.06085806194501844*sol_nodal_global[globalIdxs[2]]-0.06085806194501844*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[16] = 0.06085806194501844*sol_nodal_global[globalIdxs[19]]+0.06085806194501844*sol_nodal_global[globalIdxs[17]]-0.06085806194501844*sol_nodal_global[globalIdxs[14]]-0.06085806194501844*sol_nodal_global[globalIdxs[12]]-0.1217161238900369*sol_nodal_global[globalIdxs[11]]-0.1217161238900369*sol_nodal_global[globalIdxs[10]]+0.1217161238900369*sol_nodal_global[globalIdxs[9]]+0.1217161238900369*sol_nodal_global[globalIdxs[8]]+0.06085806194501844*sol_nodal_global[globalIdxs[7]]+0.06085806194501844*sol_nodal_global[globalIdxs[5]]-0.06085806194501844*sol_nodal_global[globalIdxs[2]]-0.06085806194501844*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[17] = 0.03513641844631532*sol_nodal_global[globalIdxs[19]]-0.07027283689263064*sol_nodal_global[globalIdxs[18]]+0.03513641844631532*sol_nodal_global[globalIdxs[17]]-0.03513641844631532*sol_nodal_global[globalIdxs[14]]+0.07027283689263064*sol_nodal_global[globalIdxs[13]]-0.03513641844631532*sol_nodal_global[globalIdxs[12]]-0.03513641844631532*sol_nodal_global[globalIdxs[7]]+0.07027283689263064*sol_nodal_global[globalIdxs[6]]-0.03513641844631532*sol_nodal_global[globalIdxs[5]]+0.03513641844631532*sol_nodal_global[globalIdxs[2]]-0.07027283689263064*sol_nodal_global[globalIdxs[1]]+0.03513641844631532*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[18] = 0.03513641844631532*sol_nodal_global[globalIdxs[19]]-0.03513641844631532*sol_nodal_global[globalIdxs[17]]-0.07027283689263064*sol_nodal_global[globalIdxs[16]]+0.07027283689263064*sol_nodal_global[globalIdxs[15]]+0.03513641844631532*sol_nodal_global[globalIdxs[14]]-0.03513641844631532*sol_nodal_global[globalIdxs[12]]-0.03513641844631532*sol_nodal_global[globalIdxs[7]]+0.03513641844631532*sol_nodal_global[globalIdxs[5]]+0.07027283689263064*sol_nodal_global[globalIdxs[4]]-0.07027283689263064*sol_nodal_global[globalIdxs[3]]-0.03513641844631532*sol_nodal_global[globalIdxs[2]]+0.03513641844631532*sol_nodal_global[globalIdxs[0]];
  sol_modal_local[19] = 0.03513641844631532*sol_nodal_global[globalIdxs[19]]-0.03513641844631532*sol_nodal_global[globalIdxs[17]]-0.03513641844631532*sol_nodal_global[globalIdxs[14]]+0.03513641844631532*sol_nodal_global[globalIdxs[12]]-0.07027283689263064*sol_nodal_global[globalIdxs[11]]+0.07027283689263064*sol_nodal_global[globalIdxs[10]]+0.07027283689263064*sol_nodal_global[globalIdxs[9]]-0.07027283689263064*sol_nodal_global[globalIdxs[8]]+0.03513641844631532*sol_nodal_global[globalIdxs[7]]-0.03513641844631532*sol_nodal_global[globalIdxs[5]]-0.03513641844631532*sol_nodal_global[globalIdxs[2]]+0.03513641844631532*sol_nodal_global[globalIdxs[0]];

}
