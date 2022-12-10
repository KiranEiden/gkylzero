#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_3x_p2_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfx_2x1v_ser_p2(const double *w, const double *dxv, 
     const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *u_il, const double *u_ic, const double *u_ir, 
     const double *T_ijl, const double *T_ijc, const double *T_ijr,
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                 Cell-center coordinates.
  // dxv[NDIM]:               Cell spacing.
  // bvarl/bvarc/bvarr:       Input magnetic field unit vector in left/center/right cells.
  // u_il/u_ic/u_ir:          Input bulk velocity (ux,uy,uz) in left/center/right cells.
  // T_ijl/T_ijc/T_ijr:       Input Temperature tensor/mass (for penalization) in left/center/right cells.
  // fl/fc/fr:                Input Distribution function in left/center/right cells.
  // out:                     Incremented distribution function in center cell.
  const double dx1 = 2.0/dxv[0]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *ul = &u_il[0]; 
  const double *uc = &u_ic[0]; 
  const double *ur = &u_ir[0]; 
  const double *bl = &bvarl[0]; 
  const double *bc = &bvarc[0]; 
  const double *br = &bvarr[0]; 
  // Get thermal velocity in direction of update for penalization vth^2 = 3.0*T_ii/m. 
  const double *vth_sql = &T_ijl[0]; 
  const double *vth_sqc = &T_ijc[0]; 
  const double *vth_sqr = &T_ijr[0]; 

  double alpha_l[20] = {0.0}; 
  double alpha_c[20] = {0.0}; 
  double alpha_r[20] = {0.0}; 
  alpha_l[0] = 1.414213562373095*bl[0]*wvpar; 
  alpha_l[1] = 1.414213562373095*bl[1]*wvpar; 
  alpha_l[2] = 1.414213562373095*bl[2]*wvpar; 
  alpha_l[3] = 0.408248290463863*bl[0]*dvpar; 
  alpha_l[4] = 1.414213562373095*bl[3]*wvpar; 
  alpha_l[5] = 0.408248290463863*bl[1]*dvpar; 
  alpha_l[6] = 0.408248290463863*bl[2]*dvpar; 
  alpha_l[7] = 1.414213562373095*bl[4]*wvpar; 
  alpha_l[8] = 1.414213562373095*bl[5]*wvpar; 
  alpha_l[10] = 0.408248290463863*bl[3]*dvpar; 
  alpha_l[11] = 1.414213562373095*bl[6]*wvpar; 
  alpha_l[12] = 1.414213562373095*bl[7]*wvpar; 
  alpha_l[13] = 0.408248290463863*bl[4]*dvpar; 
  alpha_l[14] = 0.408248290463863*bl[5]*dvpar; 
  alpha_l[17] = 0.408248290463863*bl[6]*dvpar; 
  alpha_l[18] = 0.408248290463863*bl[7]*dvpar; 

  alpha_c[0] = 1.414213562373095*bc[0]*wvpar; 
  alpha_c[1] = 1.414213562373095*bc[1]*wvpar; 
  alpha_c[2] = 1.414213562373095*bc[2]*wvpar; 
  alpha_c[3] = 0.408248290463863*bc[0]*dvpar; 
  alpha_c[4] = 1.414213562373095*bc[3]*wvpar; 
  alpha_c[5] = 0.408248290463863*bc[1]*dvpar; 
  alpha_c[6] = 0.408248290463863*bc[2]*dvpar; 
  alpha_c[7] = 1.414213562373095*bc[4]*wvpar; 
  alpha_c[8] = 1.414213562373095*bc[5]*wvpar; 
  alpha_c[10] = 0.408248290463863*bc[3]*dvpar; 
  alpha_c[11] = 1.414213562373095*bc[6]*wvpar; 
  alpha_c[12] = 1.414213562373095*bc[7]*wvpar; 
  alpha_c[13] = 0.408248290463863*bc[4]*dvpar; 
  alpha_c[14] = 0.408248290463863*bc[5]*dvpar; 
  alpha_c[17] = 0.408248290463863*bc[6]*dvpar; 
  alpha_c[18] = 0.408248290463863*bc[7]*dvpar; 

  alpha_r[0] = 1.414213562373095*br[0]*wvpar; 
  alpha_r[1] = 1.414213562373095*br[1]*wvpar; 
  alpha_r[2] = 1.414213562373095*br[2]*wvpar; 
  alpha_r[3] = 0.408248290463863*br[0]*dvpar; 
  alpha_r[4] = 1.414213562373095*br[3]*wvpar; 
  alpha_r[5] = 0.408248290463863*br[1]*dvpar; 
  alpha_r[6] = 0.408248290463863*br[2]*dvpar; 
  alpha_r[7] = 1.414213562373095*br[4]*wvpar; 
  alpha_r[8] = 1.414213562373095*br[5]*wvpar; 
  alpha_r[10] = 0.408248290463863*br[3]*dvpar; 
  alpha_r[11] = 1.414213562373095*br[6]*wvpar; 
  alpha_r[12] = 1.414213562373095*br[7]*wvpar; 
  alpha_r[13] = 0.408248290463863*br[4]*dvpar; 
  alpha_r[14] = 0.408248290463863*br[5]*dvpar; 
  alpha_r[17] = 0.408248290463863*br[6]*dvpar; 
  alpha_r[18] = 0.408248290463863*br[7]*dvpar; 

  double alphaSurf_l[8] = {0.0}; 
  alphaSurf_l[0] = 0.3458741190809163*alpha_l[7]+0.3458741190809163*alpha_c[7]+0.4975526040028326*alpha_l[1]-0.4975526040028326*alpha_c[1]+0.3535533905932737*alpha_l[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_l[1] = 0.3458741190809163*alpha_l[11]+0.3458741190809163*alpha_c[11]+0.4975526040028326*alpha_l[4]-0.4975526040028326*alpha_c[4]+0.3535533905932737*alpha_l[2]+0.3535533905932737*alpha_c[2]; 
  alphaSurf_l[2] = 0.3458741190809163*alpha_l[13]+0.3458741190809163*alpha_c[13]+0.4975526040028326*alpha_l[5]-0.4975526040028326*alpha_c[5]+0.3535533905932737*alpha_l[3]+0.3535533905932737*alpha_c[3]; 
  alphaSurf_l[3] = 0.3458741190809163*alpha_l[17]+0.3458741190809163*alpha_c[17]+0.4975526040028326*alpha_l[10]-0.4975526040028326*alpha_c[10]+0.3535533905932737*alpha_l[6]+0.3535533905932737*alpha_c[6]; 
  alphaSurf_l[4] = 0.4975526040028326*alpha_l[12]-0.4975526040028326*alpha_c[12]+0.3535533905932737*alpha_l[8]+0.3535533905932737*alpha_c[8]; 
  alphaSurf_l[6] = 0.4975526040028326*alpha_l[18]-0.4975526040028326*alpha_c[18]+0.3535533905932737*alpha_l[14]+0.3535533905932737*alpha_c[14]; 

  double alphaSurf_r[8] = {0.0}; 
  alphaSurf_r[0] = 0.3458741190809163*alpha_r[7]+0.3458741190809163*alpha_c[7]-0.4975526040028326*alpha_r[1]+0.4975526040028326*alpha_c[1]+0.3535533905932737*alpha_r[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_r[1] = 0.3458741190809163*alpha_r[11]+0.3458741190809163*alpha_c[11]-0.4975526040028326*alpha_r[4]+0.4975526040028326*alpha_c[4]+0.3535533905932737*alpha_r[2]+0.3535533905932737*alpha_c[2]; 
  alphaSurf_r[2] = 0.3458741190809163*alpha_r[13]+0.3458741190809163*alpha_c[13]-0.4975526040028326*alpha_r[5]+0.4975526040028326*alpha_c[5]+0.3535533905932737*alpha_r[3]+0.3535533905932737*alpha_c[3]; 
  alphaSurf_r[3] = 0.3458741190809163*alpha_r[17]+0.3458741190809163*alpha_c[17]-0.4975526040028326*alpha_r[10]+0.4975526040028326*alpha_c[10]+0.3535533905932737*alpha_r[6]+0.3535533905932737*alpha_c[6]; 
  alphaSurf_r[4] = (-0.4975526040028326*alpha_r[12])+0.4975526040028326*alpha_c[12]+0.3535533905932737*alpha_r[8]+0.3535533905932737*alpha_c[8]; 
  alphaSurf_r[6] = (-0.4975526040028326*alpha_r[18])+0.4975526040028326*alpha_c[18]+0.3535533905932737*alpha_r[14]+0.3535533905932737*alpha_c[14]; 

  double fUpwindQuad_l[9] = {0.0};
  double fUpwindQuad_r[9] = {0.0};
  double fUpwind_l[8] = {0.0};
  double fUpwind_r[8] = {0.0};
  double Ghat_l[8] = {0.0}; 
  double Ghat_r[8] = {0.0}; 

  if ((-0.5999999999999999*alphaSurf_l[6])+0.4472135954999579*alphaSurf_l[4]+0.9*alphaSurf_l[3]-0.6708203932499369*(alphaSurf_l[2]+alphaSurf_l[1])+0.5*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[0] = ser_3x_p2_surfx1_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_3x_p2_surfx1_eval_quad_node_0_l(fc); 
  } 
  if ((-0.5999999999999999*alphaSurf_r[6])+0.4472135954999579*alphaSurf_r[4]+0.9*alphaSurf_r[3]-0.6708203932499369*(alphaSurf_r[2]+alphaSurf_r[1])+0.5*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[0] = ser_3x_p2_surfx1_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_3x_p2_surfx1_eval_quad_node_0_l(fr); 
  } 
  if (0.4472135954999579*alphaSurf_l[4]-0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[1] = ser_3x_p2_surfx1_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_3x_p2_surfx1_eval_quad_node_1_l(fc); 
  } 
  if (0.4472135954999579*alphaSurf_r[4]-0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[1] = ser_3x_p2_surfx1_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_3x_p2_surfx1_eval_quad_node_1_l(fr); 
  } 
  if (0.5999999999999999*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[4]-0.9*alphaSurf_l[3]+0.6708203932499369*alphaSurf_l[2]-0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = ser_3x_p2_surfx1_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_3x_p2_surfx1_eval_quad_node_2_l(fc); 
  } 
  if (0.5999999999999999*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[4]-0.9*alphaSurf_r[3]+0.6708203932499369*alphaSurf_r[2]-0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[2] = ser_3x_p2_surfx1_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_3x_p2_surfx1_eval_quad_node_2_l(fr); 
  } 
  if (0.75*alphaSurf_l[6]-0.5590169943749475*alphaSurf_l[4]-0.6708203932499369*alphaSurf_l[2]+0.5*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[3] = ser_3x_p2_surfx1_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[3] = ser_3x_p2_surfx1_eval_quad_node_3_l(fc); 
  } 
  if (0.75*alphaSurf_r[6]-0.5590169943749475*alphaSurf_r[4]-0.6708203932499369*alphaSurf_r[2]+0.5*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[3] = ser_3x_p2_surfx1_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[3] = ser_3x_p2_surfx1_eval_quad_node_3_l(fr); 
  } 
  if (0.5*alphaSurf_l[0]-0.5590169943749475*alphaSurf_l[4] > 0) { 
    fUpwindQuad_l[4] = ser_3x_p2_surfx1_eval_quad_node_4_r(fl); 
  } else { 
    fUpwindQuad_l[4] = ser_3x_p2_surfx1_eval_quad_node_4_l(fc); 
  } 
  if (0.5*alphaSurf_r[0]-0.5590169943749475*alphaSurf_r[4] > 0) { 
    fUpwindQuad_r[4] = ser_3x_p2_surfx1_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_r[4] = ser_3x_p2_surfx1_eval_quad_node_4_l(fr); 
  } 
  if ((-0.75*alphaSurf_l[6])-0.5590169943749475*alphaSurf_l[4]+0.6708203932499369*alphaSurf_l[2]+0.5*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[5] = ser_3x_p2_surfx1_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_l[5] = ser_3x_p2_surfx1_eval_quad_node_5_l(fc); 
  } 
  if ((-0.75*alphaSurf_r[6])-0.5590169943749475*alphaSurf_r[4]+0.6708203932499369*alphaSurf_r[2]+0.5*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[5] = ser_3x_p2_surfx1_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_r[5] = ser_3x_p2_surfx1_eval_quad_node_5_l(fr); 
  } 
  if ((-0.5999999999999999*alphaSurf_l[6])+0.4472135954999579*alphaSurf_l[4]-0.9*alphaSurf_l[3]-0.6708203932499369*alphaSurf_l[2]+0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[6] = ser_3x_p2_surfx1_eval_quad_node_6_r(fl); 
  } else { 
    fUpwindQuad_l[6] = ser_3x_p2_surfx1_eval_quad_node_6_l(fc); 
  } 
  if ((-0.5999999999999999*alphaSurf_r[6])+0.4472135954999579*alphaSurf_r[4]-0.9*alphaSurf_r[3]-0.6708203932499369*alphaSurf_r[2]+0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[6] = ser_3x_p2_surfx1_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_r[6] = ser_3x_p2_surfx1_eval_quad_node_6_l(fr); 
  } 
  if (0.4472135954999579*alphaSurf_l[4]+0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[7] = ser_3x_p2_surfx1_eval_quad_node_7_r(fl); 
  } else { 
    fUpwindQuad_l[7] = ser_3x_p2_surfx1_eval_quad_node_7_l(fc); 
  } 
  if (0.4472135954999579*alphaSurf_r[4]+0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[7] = ser_3x_p2_surfx1_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_r[7] = ser_3x_p2_surfx1_eval_quad_node_7_l(fr); 
  } 
  if (0.5999999999999999*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[4]+0.9*alphaSurf_l[3]+0.6708203932499369*(alphaSurf_l[2]+alphaSurf_l[1])+0.5*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[8] = ser_3x_p2_surfx1_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_l[8] = ser_3x_p2_surfx1_eval_quad_node_8_l(fc); 
  } 
  if (0.5999999999999999*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[4]+0.9*alphaSurf_r[3]+0.6708203932499369*(alphaSurf_r[2]+alphaSurf_r[1])+0.5*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[8] = ser_3x_p2_surfx1_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_r[8] = ser_3x_p2_surfx1_eval_quad_node_8_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p2_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_3x_p2_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.5*alphaSurf_l[6]*fUpwind_l[6]+0.5*alphaSurf_l[4]*fUpwind_l[4]+0.5*alphaSurf_l[3]*fUpwind_l[3]+0.5*alphaSurf_l[2]*fUpwind_l[2]+0.5*alphaSurf_l[1]*fUpwind_l[1]+0.5*alphaSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.447213595499958*alphaSurf_l[3]*fUpwind_l[6]+0.447213595499958*fUpwind_l[3]*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[1]*fUpwind_l[4]+0.4472135954999579*fUpwind_l[1]*alphaSurf_l[4]+0.5*alphaSurf_l[2]*fUpwind_l[3]+0.5*fUpwind_l[2]*alphaSurf_l[3]+0.5*alphaSurf_l[0]*fUpwind_l[1]+0.5*fUpwind_l[0]*alphaSurf_l[1]; 
  Ghat_l[2] = 0.447213595499958*alphaSurf_l[3]*fUpwind_l[7]+0.5000000000000001*alphaSurf_l[4]*fUpwind_l[6]+0.5000000000000001*fUpwind_l[4]*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[2]*fUpwind_l[5]+0.5*alphaSurf_l[1]*fUpwind_l[3]+0.5*fUpwind_l[1]*alphaSurf_l[3]+0.5*alphaSurf_l[0]*fUpwind_l[2]+0.5*fUpwind_l[0]*alphaSurf_l[2]; 
  Ghat_l[3] = 0.4*alphaSurf_l[6]*fUpwind_l[7]+0.447213595499958*alphaSurf_l[2]*fUpwind_l[7]+0.447213595499958*alphaSurf_l[1]*fUpwind_l[6]+0.447213595499958*fUpwind_l[1]*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[3]*fUpwind_l[5]+0.4472135954999579*alphaSurf_l[3]*fUpwind_l[4]+0.4472135954999579*fUpwind_l[3]*alphaSurf_l[4]+0.5*alphaSurf_l[0]*fUpwind_l[3]+0.5*fUpwind_l[0]*alphaSurf_l[3]+0.5*alphaSurf_l[1]*fUpwind_l[2]+0.5*fUpwind_l[1]*alphaSurf_l[2]; 
  Ghat_l[4] = 0.31943828249997*alphaSurf_l[6]*fUpwind_l[6]+0.5000000000000001*alphaSurf_l[2]*fUpwind_l[6]+0.5000000000000001*fUpwind_l[2]*alphaSurf_l[6]+0.31943828249997*alphaSurf_l[4]*fUpwind_l[4]+0.5*alphaSurf_l[0]*fUpwind_l[4]+0.5*fUpwind_l[0]*alphaSurf_l[4]+0.4472135954999579*alphaSurf_l[3]*fUpwind_l[3]+0.4472135954999579*alphaSurf_l[1]*fUpwind_l[1]; 
  Ghat_l[5] = 0.5000000000000001*alphaSurf_l[1]*fUpwind_l[7]+0.4472135954999579*alphaSurf_l[6]*fUpwind_l[6]+0.5*alphaSurf_l[0]*fUpwind_l[5]+0.4472135954999579*alphaSurf_l[3]*fUpwind_l[3]+0.4472135954999579*alphaSurf_l[2]*fUpwind_l[2]; 
  Ghat_l[6] = 0.4*alphaSurf_l[3]*fUpwind_l[7]+0.31943828249997*alphaSurf_l[4]*fUpwind_l[6]+0.5*alphaSurf_l[0]*fUpwind_l[6]+0.4472135954999579*fUpwind_l[5]*alphaSurf_l[6]+0.31943828249997*fUpwind_l[4]*alphaSurf_l[6]+0.5*fUpwind_l[0]*alphaSurf_l[6]+0.5000000000000001*alphaSurf_l[2]*fUpwind_l[4]+0.5000000000000001*fUpwind_l[2]*alphaSurf_l[4]+0.447213595499958*alphaSurf_l[1]*fUpwind_l[3]+0.447213595499958*fUpwind_l[1]*alphaSurf_l[3]; 
  Ghat_l[7] = 0.4472135954999579*alphaSurf_l[4]*fUpwind_l[7]+0.5*alphaSurf_l[0]*fUpwind_l[7]+0.4*alphaSurf_l[3]*fUpwind_l[6]+0.4*fUpwind_l[3]*alphaSurf_l[6]+0.5000000000000001*alphaSurf_l[1]*fUpwind_l[5]+0.447213595499958*alphaSurf_l[2]*fUpwind_l[3]+0.447213595499958*fUpwind_l[2]*alphaSurf_l[3]; 

  Ghat_r[0] = 0.5*alphaSurf_r[6]*fUpwind_r[6]+0.5*alphaSurf_r[4]*fUpwind_r[4]+0.5*alphaSurf_r[3]*fUpwind_r[3]+0.5*alphaSurf_r[2]*fUpwind_r[2]+0.5*alphaSurf_r[1]*fUpwind_r[1]+0.5*alphaSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.447213595499958*alphaSurf_r[3]*fUpwind_r[6]+0.447213595499958*fUpwind_r[3]*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[1]*fUpwind_r[4]+0.4472135954999579*fUpwind_r[1]*alphaSurf_r[4]+0.5*alphaSurf_r[2]*fUpwind_r[3]+0.5*fUpwind_r[2]*alphaSurf_r[3]+0.5*alphaSurf_r[0]*fUpwind_r[1]+0.5*fUpwind_r[0]*alphaSurf_r[1]; 
  Ghat_r[2] = 0.447213595499958*alphaSurf_r[3]*fUpwind_r[7]+0.5000000000000001*alphaSurf_r[4]*fUpwind_r[6]+0.5000000000000001*fUpwind_r[4]*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[2]*fUpwind_r[5]+0.5*alphaSurf_r[1]*fUpwind_r[3]+0.5*fUpwind_r[1]*alphaSurf_r[3]+0.5*alphaSurf_r[0]*fUpwind_r[2]+0.5*fUpwind_r[0]*alphaSurf_r[2]; 
  Ghat_r[3] = 0.4*alphaSurf_r[6]*fUpwind_r[7]+0.447213595499958*alphaSurf_r[2]*fUpwind_r[7]+0.447213595499958*alphaSurf_r[1]*fUpwind_r[6]+0.447213595499958*fUpwind_r[1]*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[3]*fUpwind_r[5]+0.4472135954999579*alphaSurf_r[3]*fUpwind_r[4]+0.4472135954999579*fUpwind_r[3]*alphaSurf_r[4]+0.5*alphaSurf_r[0]*fUpwind_r[3]+0.5*fUpwind_r[0]*alphaSurf_r[3]+0.5*alphaSurf_r[1]*fUpwind_r[2]+0.5*fUpwind_r[1]*alphaSurf_r[2]; 
  Ghat_r[4] = 0.31943828249997*alphaSurf_r[6]*fUpwind_r[6]+0.5000000000000001*alphaSurf_r[2]*fUpwind_r[6]+0.5000000000000001*fUpwind_r[2]*alphaSurf_r[6]+0.31943828249997*alphaSurf_r[4]*fUpwind_r[4]+0.5*alphaSurf_r[0]*fUpwind_r[4]+0.5*fUpwind_r[0]*alphaSurf_r[4]+0.4472135954999579*alphaSurf_r[3]*fUpwind_r[3]+0.4472135954999579*alphaSurf_r[1]*fUpwind_r[1]; 
  Ghat_r[5] = 0.5000000000000001*alphaSurf_r[1]*fUpwind_r[7]+0.4472135954999579*alphaSurf_r[6]*fUpwind_r[6]+0.5*alphaSurf_r[0]*fUpwind_r[5]+0.4472135954999579*alphaSurf_r[3]*fUpwind_r[3]+0.4472135954999579*alphaSurf_r[2]*fUpwind_r[2]; 
  Ghat_r[6] = 0.4*alphaSurf_r[3]*fUpwind_r[7]+0.31943828249997*alphaSurf_r[4]*fUpwind_r[6]+0.5*alphaSurf_r[0]*fUpwind_r[6]+0.4472135954999579*fUpwind_r[5]*alphaSurf_r[6]+0.31943828249997*fUpwind_r[4]*alphaSurf_r[6]+0.5*fUpwind_r[0]*alphaSurf_r[6]+0.5000000000000001*alphaSurf_r[2]*fUpwind_r[4]+0.5000000000000001*fUpwind_r[2]*alphaSurf_r[4]+0.447213595499958*alphaSurf_r[1]*fUpwind_r[3]+0.447213595499958*fUpwind_r[1]*alphaSurf_r[3]; 
  Ghat_r[7] = 0.4472135954999579*alphaSurf_r[4]*fUpwind_r[7]+0.5*alphaSurf_r[0]*fUpwind_r[7]+0.4*alphaSurf_r[3]*fUpwind_r[6]+0.4*fUpwind_r[3]*alphaSurf_r[6]+0.5000000000000001*alphaSurf_r[1]*fUpwind_r[5]+0.447213595499958*alphaSurf_r[2]*fUpwind_r[3]+0.447213595499958*fUpwind_r[2]*alphaSurf_r[3]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx1; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx1; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx1; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx1; 
  out[4] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx1; 
  out[5] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx1; 
  out[6] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx1; 
  out[7] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dx1; 
  out[8] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx1; 
  out[9] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx1; 
  out[10] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx1; 
  out[11] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dx1; 
  out[12] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx1; 
  out[13] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*dx1; 
  out[14] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx1; 
  out[15] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx1; 
  out[16] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx1; 
  out[17] += (1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*dx1; 
  out[18] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx1; 
  out[19] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx1; 

  double alpha_u_l[20] = {0.0}; 
  double alpha_u_c[20] = {0.0}; 
  double alpha_u_r[20] = {0.0}; 
  double vth_sq_l[20] = {0.0}; 
  double vth_sq_c[20] = {0.0}; 
  double vth_sq_r[20] = {0.0}; 
  alpha_u_l[0] = 1.414213562373095*ul[0]; 
  alpha_u_l[1] = 1.414213562373095*ul[1]; 
  alpha_u_l[2] = 1.414213562373095*ul[2]; 
  alpha_u_l[4] = 1.414213562373095*ul[3]; 
  alpha_u_l[7] = 1.414213562373095*ul[4]; 
  alpha_u_l[8] = 1.414213562373095*ul[5]; 
  alpha_u_l[11] = 1.414213562373095*ul[6]; 
  alpha_u_l[12] = 1.414213562373095*ul[7]; 

  alpha_u_c[0] = 1.414213562373095*uc[0]; 
  alpha_u_c[1] = 1.414213562373095*uc[1]; 
  alpha_u_c[2] = 1.414213562373095*uc[2]; 
  alpha_u_c[4] = 1.414213562373095*uc[3]; 
  alpha_u_c[7] = 1.414213562373095*uc[4]; 
  alpha_u_c[8] = 1.414213562373095*uc[5]; 
  alpha_u_c[11] = 1.414213562373095*uc[6]; 
  alpha_u_c[12] = 1.414213562373095*uc[7]; 

  alpha_u_r[0] = 1.414213562373095*ur[0]; 
  alpha_u_r[1] = 1.414213562373095*ur[1]; 
  alpha_u_r[2] = 1.414213562373095*ur[2]; 
  alpha_u_r[4] = 1.414213562373095*ur[3]; 
  alpha_u_r[7] = 1.414213562373095*ur[4]; 
  alpha_u_r[8] = 1.414213562373095*ur[5]; 
  alpha_u_r[11] = 1.414213562373095*ur[6]; 
  alpha_u_r[12] = 1.414213562373095*ur[7]; 

  vth_sq_l[0] = 1.414213562373095*vth_sql[0]; 
  vth_sq_l[1] = 1.414213562373095*vth_sql[1]; 
  vth_sq_l[2] = 1.414213562373095*vth_sql[2]; 
  vth_sq_l[4] = 1.414213562373095*vth_sql[3]; 
  vth_sq_l[7] = 1.414213562373095*vth_sql[4]; 
  vth_sq_l[8] = 1.414213562373095*vth_sql[5]; 
  vth_sq_l[11] = 1.414213562373095*vth_sql[6]; 
  vth_sq_l[12] = 1.414213562373095*vth_sql[7]; 

  vth_sq_c[0] = 1.414213562373095*vth_sqc[0]; 
  vth_sq_c[1] = 1.414213562373095*vth_sqc[1]; 
  vth_sq_c[2] = 1.414213562373095*vth_sqc[2]; 
  vth_sq_c[4] = 1.414213562373095*vth_sqc[3]; 
  vth_sq_c[7] = 1.414213562373095*vth_sqc[4]; 
  vth_sq_c[8] = 1.414213562373095*vth_sqc[5]; 
  vth_sq_c[11] = 1.414213562373095*vth_sqc[6]; 
  vth_sq_c[12] = 1.414213562373095*vth_sqc[7]; 

  vth_sq_r[0] = 1.414213562373095*vth_sqr[0]; 
  vth_sq_r[1] = 1.414213562373095*vth_sqr[1]; 
  vth_sq_r[2] = 1.414213562373095*vth_sqr[2]; 
  vth_sq_r[4] = 1.414213562373095*vth_sqr[3]; 
  vth_sq_r[7] = 1.414213562373095*vth_sqr[4]; 
  vth_sq_r[8] = 1.414213562373095*vth_sqr[5]; 
  vth_sq_r[11] = 1.414213562373095*vth_sqr[6]; 
  vth_sq_r[12] = 1.414213562373095*vth_sqr[7]; 

  double lax_quad_l[9] = {0.0};
  double lax_quad_r[9] = {0.0};
  double lax_modal_l[8] = {0.0};
  double lax_modal_r[8] = {0.0};

  double alpha_l_r = 0.0; 
  double alpha_c_l = 0.0; 
  double alpha_c_r = 0.0; 
  double alpha_r_l = 0.0; 
  double alphaQuad_l = 0.0; 
  double alphaQuad_r = 0.0; 
  double vth_sq_l_r = 0.0; 
  double vth_sq_c_l = 0.0; 
  double vth_sq_c_r = 0.0; 
  double vth_sq_r_l = 0.0; 
  double vthQuad_l = 0.0; 
  double vthQuad_r = 0.0; 
  double max_speed_l = 0.0; 
  double max_speed_r = 0.0; 
  double f_l_r = 0.0; 
  double f_c_l = 0.0; 
  double f_c_r = 0.0; 
  double f_r_l = 0.0; 
  alpha_l_r = ser_3x_p2_surfx1_eval_quad_node_0_r(alpha_u_l); 
  alpha_c_l = ser_3x_p2_surfx1_eval_quad_node_0_l(alpha_u_c); 
  alpha_c_r = ser_3x_p2_surfx1_eval_quad_node_0_r(alpha_u_c); 
  alpha_r_l = ser_3x_p2_surfx1_eval_quad_node_0_l(alpha_u_r); 
  alphaQuad_l = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  vth_sq_l_r = ser_3x_p2_surfx1_eval_quad_node_0_r(vth_sq_l); 
  vth_sq_c_l = ser_3x_p2_surfx1_eval_quad_node_0_l(vth_sq_c); 
  vth_sq_c_r = ser_3x_p2_surfx1_eval_quad_node_0_r(vth_sq_c); 
  vth_sq_r_l = ser_3x_p2_surfx1_eval_quad_node_0_l(vth_sq_r); 
  vthQuad_l = fmax(sqrt(fabs(vth_sq_l_r)), sqrt(fabs(vth_sq_c_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sq_c_r)), sqrt(fabs(vth_sq_r_l))); 
  max_speed_l = alphaQuad_l + vthQuad_l; 
  max_speed_r = alphaQuad_r + vthQuad_r; 
  f_l_r = ser_3x_p2_surfx1_eval_quad_node_0_r(fl); 
  f_c_l = ser_3x_p2_surfx1_eval_quad_node_0_l(fc); 
  f_c_r = ser_3x_p2_surfx1_eval_quad_node_0_r(fc); 
  f_r_l = ser_3x_p2_surfx1_eval_quad_node_0_l(fr); 
  lax_quad_l[0] = 0.5*(alpha_l_r*f_l_r + alpha_c_l*f_c_l) - 0.5*max_speed_l*(f_c_l - f_l_r); 
  lax_quad_r[0] = 0.5*(alpha_c_r*f_c_r + alpha_r_l*f_r_l) - 0.5*max_speed_r*(f_r_l - f_c_r); 

  alpha_l_r = ser_3x_p2_surfx1_eval_quad_node_1_r(alpha_u_l); 
  alpha_c_l = ser_3x_p2_surfx1_eval_quad_node_1_l(alpha_u_c); 
  alpha_c_r = ser_3x_p2_surfx1_eval_quad_node_1_r(alpha_u_c); 
  alpha_r_l = ser_3x_p2_surfx1_eval_quad_node_1_l(alpha_u_r); 
  alphaQuad_l = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  vth_sq_l_r = ser_3x_p2_surfx1_eval_quad_node_1_r(vth_sq_l); 
  vth_sq_c_l = ser_3x_p2_surfx1_eval_quad_node_1_l(vth_sq_c); 
  vth_sq_c_r = ser_3x_p2_surfx1_eval_quad_node_1_r(vth_sq_c); 
  vth_sq_r_l = ser_3x_p2_surfx1_eval_quad_node_1_l(vth_sq_r); 
  vthQuad_l = fmax(sqrt(fabs(vth_sq_l_r)), sqrt(fabs(vth_sq_c_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sq_c_r)), sqrt(fabs(vth_sq_r_l))); 
  max_speed_l = alphaQuad_l + vthQuad_l; 
  max_speed_r = alphaQuad_r + vthQuad_r; 
  f_l_r = ser_3x_p2_surfx1_eval_quad_node_1_r(fl); 
  f_c_l = ser_3x_p2_surfx1_eval_quad_node_1_l(fc); 
  f_c_r = ser_3x_p2_surfx1_eval_quad_node_1_r(fc); 
  f_r_l = ser_3x_p2_surfx1_eval_quad_node_1_l(fr); 
  lax_quad_l[1] = 0.5*(alpha_l_r*f_l_r + alpha_c_l*f_c_l) - 0.5*max_speed_l*(f_c_l - f_l_r); 
  lax_quad_r[1] = 0.5*(alpha_c_r*f_c_r + alpha_r_l*f_r_l) - 0.5*max_speed_r*(f_r_l - f_c_r); 

  alpha_l_r = ser_3x_p2_surfx1_eval_quad_node_2_r(alpha_u_l); 
  alpha_c_l = ser_3x_p2_surfx1_eval_quad_node_2_l(alpha_u_c); 
  alpha_c_r = ser_3x_p2_surfx1_eval_quad_node_2_r(alpha_u_c); 
  alpha_r_l = ser_3x_p2_surfx1_eval_quad_node_2_l(alpha_u_r); 
  alphaQuad_l = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  vth_sq_l_r = ser_3x_p2_surfx1_eval_quad_node_2_r(vth_sq_l); 
  vth_sq_c_l = ser_3x_p2_surfx1_eval_quad_node_2_l(vth_sq_c); 
  vth_sq_c_r = ser_3x_p2_surfx1_eval_quad_node_2_r(vth_sq_c); 
  vth_sq_r_l = ser_3x_p2_surfx1_eval_quad_node_2_l(vth_sq_r); 
  vthQuad_l = fmax(sqrt(fabs(vth_sq_l_r)), sqrt(fabs(vth_sq_c_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sq_c_r)), sqrt(fabs(vth_sq_r_l))); 
  max_speed_l = alphaQuad_l + vthQuad_l; 
  max_speed_r = alphaQuad_r + vthQuad_r; 
  f_l_r = ser_3x_p2_surfx1_eval_quad_node_2_r(fl); 
  f_c_l = ser_3x_p2_surfx1_eval_quad_node_2_l(fc); 
  f_c_r = ser_3x_p2_surfx1_eval_quad_node_2_r(fc); 
  f_r_l = ser_3x_p2_surfx1_eval_quad_node_2_l(fr); 
  lax_quad_l[2] = 0.5*(alpha_l_r*f_l_r + alpha_c_l*f_c_l) - 0.5*max_speed_l*(f_c_l - f_l_r); 
  lax_quad_r[2] = 0.5*(alpha_c_r*f_c_r + alpha_r_l*f_r_l) - 0.5*max_speed_r*(f_r_l - f_c_r); 

  alpha_l_r = ser_3x_p2_surfx1_eval_quad_node_3_r(alpha_u_l); 
  alpha_c_l = ser_3x_p2_surfx1_eval_quad_node_3_l(alpha_u_c); 
  alpha_c_r = ser_3x_p2_surfx1_eval_quad_node_3_r(alpha_u_c); 
  alpha_r_l = ser_3x_p2_surfx1_eval_quad_node_3_l(alpha_u_r); 
  alphaQuad_l = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  vth_sq_l_r = ser_3x_p2_surfx1_eval_quad_node_3_r(vth_sq_l); 
  vth_sq_c_l = ser_3x_p2_surfx1_eval_quad_node_3_l(vth_sq_c); 
  vth_sq_c_r = ser_3x_p2_surfx1_eval_quad_node_3_r(vth_sq_c); 
  vth_sq_r_l = ser_3x_p2_surfx1_eval_quad_node_3_l(vth_sq_r); 
  vthQuad_l = fmax(sqrt(fabs(vth_sq_l_r)), sqrt(fabs(vth_sq_c_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sq_c_r)), sqrt(fabs(vth_sq_r_l))); 
  max_speed_l = alphaQuad_l + vthQuad_l; 
  max_speed_r = alphaQuad_r + vthQuad_r; 
  f_l_r = ser_3x_p2_surfx1_eval_quad_node_3_r(fl); 
  f_c_l = ser_3x_p2_surfx1_eval_quad_node_3_l(fc); 
  f_c_r = ser_3x_p2_surfx1_eval_quad_node_3_r(fc); 
  f_r_l = ser_3x_p2_surfx1_eval_quad_node_3_l(fr); 
  lax_quad_l[3] = 0.5*(alpha_l_r*f_l_r + alpha_c_l*f_c_l) - 0.5*max_speed_l*(f_c_l - f_l_r); 
  lax_quad_r[3] = 0.5*(alpha_c_r*f_c_r + alpha_r_l*f_r_l) - 0.5*max_speed_r*(f_r_l - f_c_r); 

  alpha_l_r = ser_3x_p2_surfx1_eval_quad_node_4_r(alpha_u_l); 
  alpha_c_l = ser_3x_p2_surfx1_eval_quad_node_4_l(alpha_u_c); 
  alpha_c_r = ser_3x_p2_surfx1_eval_quad_node_4_r(alpha_u_c); 
  alpha_r_l = ser_3x_p2_surfx1_eval_quad_node_4_l(alpha_u_r); 
  alphaQuad_l = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  vth_sq_l_r = ser_3x_p2_surfx1_eval_quad_node_4_r(vth_sq_l); 
  vth_sq_c_l = ser_3x_p2_surfx1_eval_quad_node_4_l(vth_sq_c); 
  vth_sq_c_r = ser_3x_p2_surfx1_eval_quad_node_4_r(vth_sq_c); 
  vth_sq_r_l = ser_3x_p2_surfx1_eval_quad_node_4_l(vth_sq_r); 
  vthQuad_l = fmax(sqrt(fabs(vth_sq_l_r)), sqrt(fabs(vth_sq_c_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sq_c_r)), sqrt(fabs(vth_sq_r_l))); 
  max_speed_l = alphaQuad_l + vthQuad_l; 
  max_speed_r = alphaQuad_r + vthQuad_r; 
  f_l_r = ser_3x_p2_surfx1_eval_quad_node_4_r(fl); 
  f_c_l = ser_3x_p2_surfx1_eval_quad_node_4_l(fc); 
  f_c_r = ser_3x_p2_surfx1_eval_quad_node_4_r(fc); 
  f_r_l = ser_3x_p2_surfx1_eval_quad_node_4_l(fr); 
  lax_quad_l[4] = 0.5*(alpha_l_r*f_l_r + alpha_c_l*f_c_l) - 0.5*max_speed_l*(f_c_l - f_l_r); 
  lax_quad_r[4] = 0.5*(alpha_c_r*f_c_r + alpha_r_l*f_r_l) - 0.5*max_speed_r*(f_r_l - f_c_r); 

  alpha_l_r = ser_3x_p2_surfx1_eval_quad_node_5_r(alpha_u_l); 
  alpha_c_l = ser_3x_p2_surfx1_eval_quad_node_5_l(alpha_u_c); 
  alpha_c_r = ser_3x_p2_surfx1_eval_quad_node_5_r(alpha_u_c); 
  alpha_r_l = ser_3x_p2_surfx1_eval_quad_node_5_l(alpha_u_r); 
  alphaQuad_l = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  vth_sq_l_r = ser_3x_p2_surfx1_eval_quad_node_5_r(vth_sq_l); 
  vth_sq_c_l = ser_3x_p2_surfx1_eval_quad_node_5_l(vth_sq_c); 
  vth_sq_c_r = ser_3x_p2_surfx1_eval_quad_node_5_r(vth_sq_c); 
  vth_sq_r_l = ser_3x_p2_surfx1_eval_quad_node_5_l(vth_sq_r); 
  vthQuad_l = fmax(sqrt(fabs(vth_sq_l_r)), sqrt(fabs(vth_sq_c_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sq_c_r)), sqrt(fabs(vth_sq_r_l))); 
  max_speed_l = alphaQuad_l + vthQuad_l; 
  max_speed_r = alphaQuad_r + vthQuad_r; 
  f_l_r = ser_3x_p2_surfx1_eval_quad_node_5_r(fl); 
  f_c_l = ser_3x_p2_surfx1_eval_quad_node_5_l(fc); 
  f_c_r = ser_3x_p2_surfx1_eval_quad_node_5_r(fc); 
  f_r_l = ser_3x_p2_surfx1_eval_quad_node_5_l(fr); 
  lax_quad_l[5] = 0.5*(alpha_l_r*f_l_r + alpha_c_l*f_c_l) - 0.5*max_speed_l*(f_c_l - f_l_r); 
  lax_quad_r[5] = 0.5*(alpha_c_r*f_c_r + alpha_r_l*f_r_l) - 0.5*max_speed_r*(f_r_l - f_c_r); 

  alpha_l_r = ser_3x_p2_surfx1_eval_quad_node_6_r(alpha_u_l); 
  alpha_c_l = ser_3x_p2_surfx1_eval_quad_node_6_l(alpha_u_c); 
  alpha_c_r = ser_3x_p2_surfx1_eval_quad_node_6_r(alpha_u_c); 
  alpha_r_l = ser_3x_p2_surfx1_eval_quad_node_6_l(alpha_u_r); 
  alphaQuad_l = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  vth_sq_l_r = ser_3x_p2_surfx1_eval_quad_node_6_r(vth_sq_l); 
  vth_sq_c_l = ser_3x_p2_surfx1_eval_quad_node_6_l(vth_sq_c); 
  vth_sq_c_r = ser_3x_p2_surfx1_eval_quad_node_6_r(vth_sq_c); 
  vth_sq_r_l = ser_3x_p2_surfx1_eval_quad_node_6_l(vth_sq_r); 
  vthQuad_l = fmax(sqrt(fabs(vth_sq_l_r)), sqrt(fabs(vth_sq_c_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sq_c_r)), sqrt(fabs(vth_sq_r_l))); 
  max_speed_l = alphaQuad_l + vthQuad_l; 
  max_speed_r = alphaQuad_r + vthQuad_r; 
  f_l_r = ser_3x_p2_surfx1_eval_quad_node_6_r(fl); 
  f_c_l = ser_3x_p2_surfx1_eval_quad_node_6_l(fc); 
  f_c_r = ser_3x_p2_surfx1_eval_quad_node_6_r(fc); 
  f_r_l = ser_3x_p2_surfx1_eval_quad_node_6_l(fr); 
  lax_quad_l[6] = 0.5*(alpha_l_r*f_l_r + alpha_c_l*f_c_l) - 0.5*max_speed_l*(f_c_l - f_l_r); 
  lax_quad_r[6] = 0.5*(alpha_c_r*f_c_r + alpha_r_l*f_r_l) - 0.5*max_speed_r*(f_r_l - f_c_r); 

  alpha_l_r = ser_3x_p2_surfx1_eval_quad_node_7_r(alpha_u_l); 
  alpha_c_l = ser_3x_p2_surfx1_eval_quad_node_7_l(alpha_u_c); 
  alpha_c_r = ser_3x_p2_surfx1_eval_quad_node_7_r(alpha_u_c); 
  alpha_r_l = ser_3x_p2_surfx1_eval_quad_node_7_l(alpha_u_r); 
  alphaQuad_l = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  vth_sq_l_r = ser_3x_p2_surfx1_eval_quad_node_7_r(vth_sq_l); 
  vth_sq_c_l = ser_3x_p2_surfx1_eval_quad_node_7_l(vth_sq_c); 
  vth_sq_c_r = ser_3x_p2_surfx1_eval_quad_node_7_r(vth_sq_c); 
  vth_sq_r_l = ser_3x_p2_surfx1_eval_quad_node_7_l(vth_sq_r); 
  vthQuad_l = fmax(sqrt(fabs(vth_sq_l_r)), sqrt(fabs(vth_sq_c_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sq_c_r)), sqrt(fabs(vth_sq_r_l))); 
  max_speed_l = alphaQuad_l + vthQuad_l; 
  max_speed_r = alphaQuad_r + vthQuad_r; 
  f_l_r = ser_3x_p2_surfx1_eval_quad_node_7_r(fl); 
  f_c_l = ser_3x_p2_surfx1_eval_quad_node_7_l(fc); 
  f_c_r = ser_3x_p2_surfx1_eval_quad_node_7_r(fc); 
  f_r_l = ser_3x_p2_surfx1_eval_quad_node_7_l(fr); 
  lax_quad_l[7] = 0.5*(alpha_l_r*f_l_r + alpha_c_l*f_c_l) - 0.5*max_speed_l*(f_c_l - f_l_r); 
  lax_quad_r[7] = 0.5*(alpha_c_r*f_c_r + alpha_r_l*f_r_l) - 0.5*max_speed_r*(f_r_l - f_c_r); 

  alpha_l_r = ser_3x_p2_surfx1_eval_quad_node_8_r(alpha_u_l); 
  alpha_c_l = ser_3x_p2_surfx1_eval_quad_node_8_l(alpha_u_c); 
  alpha_c_r = ser_3x_p2_surfx1_eval_quad_node_8_r(alpha_u_c); 
  alpha_r_l = ser_3x_p2_surfx1_eval_quad_node_8_l(alpha_u_r); 
  alphaQuad_l = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  vth_sq_l_r = ser_3x_p2_surfx1_eval_quad_node_8_r(vth_sq_l); 
  vth_sq_c_l = ser_3x_p2_surfx1_eval_quad_node_8_l(vth_sq_c); 
  vth_sq_c_r = ser_3x_p2_surfx1_eval_quad_node_8_r(vth_sq_c); 
  vth_sq_r_l = ser_3x_p2_surfx1_eval_quad_node_8_l(vth_sq_r); 
  vthQuad_l = fmax(sqrt(fabs(vth_sq_l_r)), sqrt(fabs(vth_sq_c_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sq_c_r)), sqrt(fabs(vth_sq_r_l))); 
  max_speed_l = alphaQuad_l + vthQuad_l; 
  max_speed_r = alphaQuad_r + vthQuad_r; 
  f_l_r = ser_3x_p2_surfx1_eval_quad_node_8_r(fl); 
  f_c_l = ser_3x_p2_surfx1_eval_quad_node_8_l(fc); 
  f_c_r = ser_3x_p2_surfx1_eval_quad_node_8_r(fc); 
  f_r_l = ser_3x_p2_surfx1_eval_quad_node_8_l(fr); 
  lax_quad_l[8] = 0.5*(alpha_l_r*f_l_r + alpha_c_l*f_c_l) - 0.5*max_speed_l*(f_c_l - f_l_r); 
  lax_quad_r[8] = 0.5*(alpha_c_r*f_c_r + alpha_r_l*f_r_l) - 0.5*max_speed_r*(f_r_l - f_c_r); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p2_upwind_quad_to_modal(lax_quad_l, lax_modal_l); 
  ser_3x_p2_upwind_quad_to_modal(lax_quad_r, lax_modal_r); 

  out[0] += (0.7071067811865475*lax_modal_l[0]-0.7071067811865475*lax_modal_r[0])*dx1; 
  out[1] += -1.224744871391589*(lax_modal_r[0]+lax_modal_l[0])*dx1; 
  out[2] += (0.7071067811865475*lax_modal_l[1]-0.7071067811865475*lax_modal_r[1])*dx1; 
  out[3] += (0.7071067811865475*lax_modal_l[2]-0.7071067811865475*lax_modal_r[2])*dx1; 
  out[4] += -1.224744871391589*(lax_modal_r[1]+lax_modal_l[1])*dx1; 
  out[5] += -1.224744871391589*(lax_modal_r[2]+lax_modal_l[2])*dx1; 
  out[6] += (0.7071067811865475*lax_modal_l[3]-0.7071067811865475*lax_modal_r[3])*dx1; 
  out[7] += (1.58113883008419*lax_modal_l[0]-1.58113883008419*lax_modal_r[0])*dx1; 
  out[8] += (0.7071067811865475*lax_modal_l[4]-0.7071067811865475*lax_modal_r[4])*dx1; 
  out[9] += (0.7071067811865475*lax_modal_l[5]-0.7071067811865475*lax_modal_r[5])*dx1; 
  out[10] += -1.224744871391589*(lax_modal_r[3]+lax_modal_l[3])*dx1; 
  out[11] += (1.58113883008419*lax_modal_l[1]-1.58113883008419*lax_modal_r[1])*dx1; 
  out[12] += -1.224744871391589*(lax_modal_r[4]+lax_modal_l[4])*dx1; 
  out[13] += (1.58113883008419*lax_modal_l[2]-1.58113883008419*lax_modal_r[2])*dx1; 
  out[14] += (0.7071067811865475*lax_modal_l[6]-0.7071067811865475*lax_modal_r[6])*dx1; 
  out[15] += -1.224744871391589*(lax_modal_r[5]+lax_modal_l[5])*dx1; 
  out[16] += (0.7071067811865475*lax_modal_l[7]-0.7071067811865475*lax_modal_r[7])*dx1; 
  out[17] += (1.58113883008419*lax_modal_l[3]-1.58113883008419*lax_modal_r[3])*dx1; 
  out[18] += -1.224744871391589*(lax_modal_r[6]+lax_modal_l[6])*dx1; 
  out[19] += -1.224744871391589*(lax_modal_r[7]+lax_modal_l[7])*dx1; 


} 
