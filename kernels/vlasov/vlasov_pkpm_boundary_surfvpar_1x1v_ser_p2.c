#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_boundary_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, 
     const double *u_i, const double *div_p, const double *bvar, const double *rho_inv_b, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_i:      bulk flow velocity (ux, uy, uz).
  // div_p:     divergence of the pressure tensor.
  // bvar:      magnetic field unit vector (nine components; first three components, b_i, other six components, b_i b_j.) 
  // rho_inv_b: b_i/rho (for pressure force 1/rho * b . div(P)).
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:       Incremented distribution function in center cell.
  const double dx0 = 2.0/dxv[0]; 
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *ux = &u_i[0]; 
  const double *uy = &u_i[3]; 
  const double *uz = &u_i[6]; 

  const double *div_p_x = &div_p[0]; 
  const double *div_p_y = &div_p[3]; 
  const double *div_p_z = &div_p[6]; 

  const double *bx = &bvar[0]; 
  const double *by = &bvar[3]; 
  const double *bz = &bvar[6]; 
  const double *bxbx = &bvar[9]; 
  const double *bxby = &bvar[12]; 
  const double *bxbz = &bvar[15]; 
  const double *byby = &bvar[18]; 
  const double *bybz = &bvar[21]; 
  const double *bzbz = &bvar[24]; 

  const double *rho_inv_bx = &rho_inv_b[0]; 
  const double *rho_inv_by = &rho_inv_b[3]; 
  const double *rho_inv_bz = &rho_inv_b[6]; 

  double alphaSurf[3] = {0.0}; 
  double fUpwindQuad[3] = {0.0};
  double fUpwind[3] = {0.0};;
  double Ghat[3] = {0.0}; 

  if (edge == -1) { 

  alphaSurf[0] = ((-2.738612787525831*bxbz[1]*uz[2])-2.738612787525831*bxby[1]*uy[2]-2.738612787525831*bxbx[1]*ux[2]-1.224744871391589*bxbz[0]*uz[1]-1.224744871391589*bxby[0]*uy[1]-1.224744871391589*bxbx[0]*ux[1])*dx0*wvpar+((-1.369306393762915*bxbz[1]*uz[2])-1.369306393762915*bxby[1]*uy[2]-1.369306393762915*bxbx[1]*ux[2]-0.6123724356957944*bxbz[0]*uz[1]-0.6123724356957944*bxby[0]*uy[1]-0.6123724356957944*bxbx[0]*ux[1])*dvpar*dx0+0.7071067811865475*div_p_z[2]*rho_inv_bz[2]+0.7071067811865475*div_p_y[2]*rho_inv_by[2]+0.7071067811865475*div_p_x[2]*rho_inv_bx[2]+0.7071067811865475*div_p_z[1]*rho_inv_bz[1]+0.7071067811865475*div_p_y[1]*rho_inv_by[1]+0.7071067811865475*div_p_x[1]*rho_inv_bx[1]+0.7071067811865475*div_p_z[0]*rho_inv_bz[0]+0.7071067811865475*div_p_y[0]*rho_inv_by[0]+0.7071067811865475*div_p_x[0]*rho_inv_bx[0]; 
  alphaSurf[1] = ((-2.449489742783178*bxbz[2]*uz[2])-2.738612787525831*bxbz[0]*uz[2]-2.449489742783178*bxby[2]*uy[2]-2.738612787525831*bxby[0]*uy[2]-2.449489742783178*bxbx[2]*ux[2]-2.738612787525831*bxbx[0]*ux[2]-1.224744871391589*bxbz[1]*uz[1]-1.224744871391589*bxby[1]*uy[1]-1.224744871391589*bxbx[1]*ux[1])*dx0*wvpar+((-1.224744871391589*bxbz[2]*uz[2])-1.369306393762915*bxbz[0]*uz[2]-1.224744871391589*bxby[2]*uy[2]-1.369306393762915*bxby[0]*uy[2]-1.224744871391589*bxbx[2]*ux[2]-1.369306393762915*bxbx[0]*ux[2]-0.6123724356957944*bxbz[1]*uz[1]-0.6123724356957944*bxby[1]*uy[1]-0.6123724356957944*bxbx[1]*ux[1])*dvpar*dx0+0.6324555320336759*div_p_z[1]*rho_inv_bz[2]+0.6324555320336759*div_p_y[1]*rho_inv_by[2]+0.6324555320336759*div_p_x[1]*rho_inv_bx[2]+0.6324555320336759*rho_inv_bz[1]*div_p_z[2]+0.6324555320336759*rho_inv_by[1]*div_p_y[2]+0.6324555320336759*rho_inv_bx[1]*div_p_x[2]+0.7071067811865475*div_p_z[0]*rho_inv_bz[1]+0.7071067811865475*div_p_y[0]*rho_inv_by[1]+0.7071067811865475*div_p_x[0]*rho_inv_bx[1]+0.7071067811865475*rho_inv_bz[0]*div_p_z[1]+0.7071067811865475*rho_inv_by[0]*div_p_y[1]+0.7071067811865475*rho_inv_bx[0]*div_p_x[1]; 
  alphaSurf[2] = ((-2.449489742783178*bxbz[1]*uz[2])-2.449489742783178*bxby[1]*uy[2]-2.449489742783178*bxbx[1]*ux[2]-1.224744871391589*uz[1]*bxbz[2]-1.224744871391589*uy[1]*bxby[2]-1.224744871391589*ux[1]*bxbx[2])*dx0*wvpar+((-1.224744871391589*bxbz[1]*uz[2])-1.224744871391589*bxby[1]*uy[2]-1.224744871391589*bxbx[1]*ux[2]-0.6123724356957944*uz[1]*bxbz[2]-0.6123724356957944*uy[1]*bxby[2]-0.6123724356957944*ux[1]*bxbx[2])*dvpar*dx0+0.4517539514526256*div_p_z[2]*rho_inv_bz[2]+0.7071067811865475*div_p_z[0]*rho_inv_bz[2]+0.4517539514526256*div_p_y[2]*rho_inv_by[2]+0.7071067811865475*div_p_y[0]*rho_inv_by[2]+0.4517539514526256*div_p_x[2]*rho_inv_bx[2]+0.7071067811865475*div_p_x[0]*rho_inv_bx[2]+0.7071067811865475*rho_inv_bz[0]*div_p_z[2]+0.7071067811865475*rho_inv_by[0]*div_p_y[2]+0.7071067811865475*rho_inv_bx[0]*div_p_x[2]+0.6324555320336759*div_p_z[1]*rho_inv_bz[1]+0.6324555320336759*div_p_y[1]*rho_inv_by[1]+0.6324555320336759*div_p_x[1]*rho_inv_bx[1]; 

  if (0.6324555320336759*alphaSurf[2]-0.9486832980505137*alphaSurf[1]+0.7071067811865475*alphaSurf[0] > 0) { 
    fUpwindQuad[0] = ser_2x_p2_surfx2_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_2x_p2_surfx2_eval_quad_node_0_l(fEdge); 
  } 
  if (0.7071067811865475*alphaSurf[0]-0.7905694150420947*alphaSurf[2] > 0) { 
    fUpwindQuad[1] = ser_2x_p2_surfx2_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_2x_p2_surfx2_eval_quad_node_1_l(fEdge); 
  } 
  if (0.6324555320336759*alphaSurf[2]+0.9486832980505137*alphaSurf[1]+0.7071067811865475*alphaSurf[0] > 0) { 
    fUpwindQuad[2] = ser_2x_p2_surfx2_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_2x_p2_surfx2_eval_quad_node_2_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*alphaSurf[2]*fUpwind[2]+0.7071067811865475*alphaSurf[1]*fUpwind[1]+0.7071067811865475*alphaSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.6324555320336759*alphaSurf[1]*fUpwind[2]+0.6324555320336759*fUpwind[1]*alphaSurf[2]+0.7071067811865475*alphaSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaSurf[1]; 
  Ghat[2] = 0.4517539514526256*alphaSurf[2]*fUpwind[2]+0.7071067811865475*alphaSurf[0]*fUpwind[2]+0.7071067811865475*fUpwind[0]*alphaSurf[2]+0.6324555320336759*alphaSurf[1]*fUpwind[1]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += -0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += -1.224744871391589*Ghat[0]*dv1par; 
  out[3] += -1.224744871391589*Ghat[1]*dv1par; 
  out[4] += -0.7071067811865475*Ghat[2]*dv1par; 
  out[5] += -1.58113883008419*Ghat[0]*dv1par; 
  out[6] += -1.224744871391589*Ghat[2]*dv1par; 
  out[7] += -1.58113883008419*Ghat[1]*dv1par; 

  } else { 

  alphaSurf[0] = ((-2.738612787525831*bxbz[1]*uz[2])-2.738612787525831*bxby[1]*uy[2]-2.738612787525831*bxbx[1]*ux[2]-1.224744871391589*bxbz[0]*uz[1]-1.224744871391589*bxby[0]*uy[1]-1.224744871391589*bxbx[0]*ux[1])*dx0*wvpar+(1.369306393762915*bxbz[1]*uz[2]+1.369306393762915*bxby[1]*uy[2]+1.369306393762915*bxbx[1]*ux[2]+0.6123724356957944*bxbz[0]*uz[1]+0.6123724356957944*bxby[0]*uy[1]+0.6123724356957944*bxbx[0]*ux[1])*dvpar*dx0+0.7071067811865475*div_p_z[2]*rho_inv_bz[2]+0.7071067811865475*div_p_y[2]*rho_inv_by[2]+0.7071067811865475*div_p_x[2]*rho_inv_bx[2]+0.7071067811865475*div_p_z[1]*rho_inv_bz[1]+0.7071067811865475*div_p_y[1]*rho_inv_by[1]+0.7071067811865475*div_p_x[1]*rho_inv_bx[1]+0.7071067811865475*div_p_z[0]*rho_inv_bz[0]+0.7071067811865475*div_p_y[0]*rho_inv_by[0]+0.7071067811865475*div_p_x[0]*rho_inv_bx[0]; 
  alphaSurf[1] = ((-2.449489742783178*bxbz[2]*uz[2])-2.738612787525831*bxbz[0]*uz[2]-2.449489742783178*bxby[2]*uy[2]-2.738612787525831*bxby[0]*uy[2]-2.449489742783178*bxbx[2]*ux[2]-2.738612787525831*bxbx[0]*ux[2]-1.224744871391589*bxbz[1]*uz[1]-1.224744871391589*bxby[1]*uy[1]-1.224744871391589*bxbx[1]*ux[1])*dx0*wvpar+(1.224744871391589*bxbz[2]*uz[2]+1.369306393762915*bxbz[0]*uz[2]+1.224744871391589*bxby[2]*uy[2]+1.369306393762915*bxby[0]*uy[2]+1.224744871391589*bxbx[2]*ux[2]+1.369306393762915*bxbx[0]*ux[2]+0.6123724356957944*bxbz[1]*uz[1]+0.6123724356957944*bxby[1]*uy[1]+0.6123724356957944*bxbx[1]*ux[1])*dvpar*dx0+0.6324555320336759*div_p_z[1]*rho_inv_bz[2]+0.6324555320336759*div_p_y[1]*rho_inv_by[2]+0.6324555320336759*div_p_x[1]*rho_inv_bx[2]+0.6324555320336759*rho_inv_bz[1]*div_p_z[2]+0.6324555320336759*rho_inv_by[1]*div_p_y[2]+0.6324555320336759*rho_inv_bx[1]*div_p_x[2]+0.7071067811865475*div_p_z[0]*rho_inv_bz[1]+0.7071067811865475*div_p_y[0]*rho_inv_by[1]+0.7071067811865475*div_p_x[0]*rho_inv_bx[1]+0.7071067811865475*rho_inv_bz[0]*div_p_z[1]+0.7071067811865475*rho_inv_by[0]*div_p_y[1]+0.7071067811865475*rho_inv_bx[0]*div_p_x[1]; 
  alphaSurf[2] = ((-2.449489742783178*bxbz[1]*uz[2])-2.449489742783178*bxby[1]*uy[2]-2.449489742783178*bxbx[1]*ux[2]-1.224744871391589*uz[1]*bxbz[2]-1.224744871391589*uy[1]*bxby[2]-1.224744871391589*ux[1]*bxbx[2])*dx0*wvpar+(1.224744871391589*bxbz[1]*uz[2]+1.224744871391589*bxby[1]*uy[2]+1.224744871391589*bxbx[1]*ux[2]+0.6123724356957944*uz[1]*bxbz[2]+0.6123724356957944*uy[1]*bxby[2]+0.6123724356957944*ux[1]*bxbx[2])*dvpar*dx0+0.4517539514526256*div_p_z[2]*rho_inv_bz[2]+0.7071067811865475*div_p_z[0]*rho_inv_bz[2]+0.4517539514526256*div_p_y[2]*rho_inv_by[2]+0.7071067811865475*div_p_y[0]*rho_inv_by[2]+0.4517539514526256*div_p_x[2]*rho_inv_bx[2]+0.7071067811865475*div_p_x[0]*rho_inv_bx[2]+0.7071067811865475*rho_inv_bz[0]*div_p_z[2]+0.7071067811865475*rho_inv_by[0]*div_p_y[2]+0.7071067811865475*rho_inv_bx[0]*div_p_x[2]+0.6324555320336759*div_p_z[1]*rho_inv_bz[1]+0.6324555320336759*div_p_y[1]*rho_inv_by[1]+0.6324555320336759*div_p_x[1]*rho_inv_bx[1]; 

  if (0.6324555320336759*alphaSurf[2]-0.9486832980505137*alphaSurf[1]+0.7071067811865475*alphaSurf[0] > 0) { 
    fUpwindQuad[0] = ser_2x_p2_surfx2_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_2x_p2_surfx2_eval_quad_node_0_l(fSkin); 
  } 
  if (0.7071067811865475*alphaSurf[0]-0.7905694150420947*alphaSurf[2] > 0) { 
    fUpwindQuad[1] = ser_2x_p2_surfx2_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_2x_p2_surfx2_eval_quad_node_1_l(fSkin); 
  } 
  if (0.6324555320336759*alphaSurf[2]+0.9486832980505137*alphaSurf[1]+0.7071067811865475*alphaSurf[0] > 0) { 
    fUpwindQuad[2] = ser_2x_p2_surfx2_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_2x_p2_surfx2_eval_quad_node_2_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*alphaSurf[2]*fUpwind[2]+0.7071067811865475*alphaSurf[1]*fUpwind[1]+0.7071067811865475*alphaSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.6324555320336759*alphaSurf[1]*fUpwind[2]+0.6324555320336759*fUpwind[1]*alphaSurf[2]+0.7071067811865475*alphaSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaSurf[1]; 
  Ghat[2] = 0.4517539514526256*alphaSurf[2]*fUpwind[2]+0.7071067811865475*alphaSurf[0]*fUpwind[2]+0.7071067811865475*fUpwind[0]*alphaSurf[2]+0.6324555320336759*alphaSurf[1]*fUpwind[1]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += 0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += -1.224744871391589*Ghat[0]*dv1par; 
  out[3] += -1.224744871391589*Ghat[1]*dv1par; 
  out[4] += 0.7071067811865475*Ghat[2]*dv1par; 
  out[5] += 1.58113883008419*Ghat[0]*dv1par; 
  out[6] += -1.224744871391589*Ghat[2]*dv1par; 
  out[7] += 1.58113883008419*Ghat[1]*dv1par; 

  } 
} 
