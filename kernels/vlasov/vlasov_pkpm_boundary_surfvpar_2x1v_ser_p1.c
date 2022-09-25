#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_2x1v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_2x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_boundary_surfvpar_2x1v_ser_p1(const double *w, const double *dxv, 
     const double *u_i, const double *p_ij, const double *bvar, const double *rho_inv_b, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u_i:      bulk flow velocity (ux, uy, uz).
  // p_ij:     pressure tensor (P_xx, P_xy, P_xz, P_yy, P_yz, P_zz).
  // bvar:      magnetic field unit vector (nine components; first three components, b_i, other six components, b_i b_j.) 
  // rho_inv_b: b_i/rho (for pressure force 1/rho * b . div(P)).
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:       Incremented distribution function in center cell.
  const double dx0 = 2.0/dxv[0]; 
  const double dx1 = 2.0/dxv[1]; 
  const double dv1par = 2.0/dxv[2]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *ux = &u_i[0]; 
  const double *uy = &u_i[4]; 
  const double *uz = &u_i[8]; 

  const double *Pxx = &p_ij[0]; 
  const double *Pxy = &p_ij[4]; 
  const double *Pxz = &p_ij[8]; 
  const double *Pyy = &p_ij[12]; 
  const double *Pyz = &p_ij[16]; 
  const double *Pzz = &p_ij[20]; 

  const double *bx = &bvar[0]; 
  const double *by = &bvar[4]; 
  const double *bz = &bvar[8]; 
  const double *bxbx = &bvar[12]; 
  const double *bxby = &bvar[16]; 
  const double *bxbz = &bvar[20]; 
  const double *byby = &bvar[24]; 
  const double *bybz = &bvar[28]; 
  const double *bzbz = &bvar[32]; 

  const double *rho_inv_bx = &rho_inv_b[0]; 
  const double *rho_inv_by = &rho_inv_b[4]; 
  const double *rho_inv_bz = &rho_inv_b[8]; 

  double alphaSurf[4] = {0.0}; 
  double fUpwindQuad[4] = {0.0};
  double fUpwind[4] = {0.0};;
  double Ghat[4] = {0.0}; 

  if (edge == -1) { 

  alphaSurf[0] = ((-0.8660254037844386*bybz[1]*uz[3])-0.8660254037844386*byby[1]*uy[3]-0.8660254037844386*bxby[1]*ux[3]-0.8660254037844386*bybz[0]*uz[2]-0.8660254037844386*byby[0]*uy[2]-0.8660254037844386*bxby[0]*ux[2])*dx1*wvpar+((-0.8660254037844386*bxbz[2]*uz[3])-0.8660254037844386*bxby[2]*uy[3]-0.8660254037844386*bxbx[2]*ux[3]-0.8660254037844386*bxbz[0]*uz[1]-0.8660254037844386*bxby[0]*uy[1]-0.8660254037844386*bxbx[0]*ux[1])*dx0*wvpar+((-0.4330127018922193*bybz[1]*uz[3])-0.4330127018922193*byby[1]*uy[3]-0.4330127018922193*bxby[1]*ux[3]-0.4330127018922193*bybz[0]*uz[2]-0.4330127018922193*byby[0]*uy[2]-0.4330127018922193*bxby[0]*ux[2])*dvpar*dx1+(0.8660254037844386*rho_inv_bz[1]*Pyz[3]+0.8660254037844386*rho_inv_by[1]*Pyy[3]+0.8660254037844386*rho_inv_bx[1]*Pxy[3]+0.8660254037844386*rho_inv_bz[0]*Pyz[2]+0.8660254037844386*rho_inv_by[0]*Pyy[2]+0.8660254037844386*rho_inv_bx[0]*Pxy[2])*dx1+((-0.4330127018922193*bxbz[2]*uz[3])-0.4330127018922193*bxby[2]*uy[3]-0.4330127018922193*bxbx[2]*ux[3]-0.4330127018922193*bxbz[0]*uz[1]-0.4330127018922193*bxby[0]*uy[1]-0.4330127018922193*bxbx[0]*ux[1])*dvpar*dx0+(0.8660254037844386*rho_inv_bz[2]*Pxz[3]+0.8660254037844386*rho_inv_by[2]*Pxy[3]+0.8660254037844386*rho_inv_bx[2]*Pxx[3]+0.8660254037844386*rho_inv_bz[0]*Pxz[1]+0.8660254037844386*rho_inv_by[0]*Pxy[1]+0.8660254037844386*rho_inv_bx[0]*Pxx[1])*dx0; 
  alphaSurf[1] = ((-0.8660254037844386*bybz[0]*uz[3])-0.8660254037844386*byby[0]*uy[3]-0.8660254037844386*bxby[0]*ux[3]-0.8660254037844386*bybz[1]*uz[2]-0.8660254037844386*byby[1]*uy[2]-0.8660254037844386*bxby[1]*ux[2])*dx1*wvpar+((-0.8660254037844386*bxbz[3]*uz[3])-0.8660254037844386*bxby[3]*uy[3]-0.8660254037844386*bxbx[3]*ux[3]-0.8660254037844386*bxbz[1]*uz[1]-0.8660254037844386*bxby[1]*uy[1]-0.8660254037844386*bxbx[1]*ux[1])*dx0*wvpar+((-0.4330127018922193*bybz[0]*uz[3])-0.4330127018922193*byby[0]*uy[3]-0.4330127018922193*bxby[0]*ux[3]-0.4330127018922193*bybz[1]*uz[2]-0.4330127018922193*byby[1]*uy[2]-0.4330127018922193*bxby[1]*ux[2])*dvpar*dx1+(0.8660254037844386*rho_inv_bz[0]*Pyz[3]+0.8660254037844386*rho_inv_by[0]*Pyy[3]+0.8660254037844386*rho_inv_bx[0]*Pxy[3]+0.8660254037844386*rho_inv_bz[1]*Pyz[2]+0.8660254037844386*rho_inv_by[1]*Pyy[2]+0.8660254037844386*rho_inv_bx[1]*Pxy[2])*dx1+((-0.4330127018922193*bxbz[3]*uz[3])-0.4330127018922193*bxby[3]*uy[3]-0.4330127018922193*bxbx[3]*ux[3]-0.4330127018922193*bxbz[1]*uz[1]-0.4330127018922193*bxby[1]*uy[1]-0.4330127018922193*bxbx[1]*ux[1])*dvpar*dx0+(0.8660254037844386*Pxz[3]*rho_inv_bz[3]+0.8660254037844386*Pxy[3]*rho_inv_by[3]+0.8660254037844386*Pxx[3]*rho_inv_bx[3]+0.8660254037844386*Pxz[1]*rho_inv_bz[1]+0.8660254037844386*Pxy[1]*rho_inv_by[1]+0.8660254037844386*Pxx[1]*rho_inv_bx[1])*dx0; 
  alphaSurf[2] = ((-0.8660254037844386*bybz[3]*uz[3])-0.8660254037844386*byby[3]*uy[3]-0.8660254037844386*bxby[3]*ux[3]-0.8660254037844386*bybz[2]*uz[2]-0.8660254037844386*byby[2]*uy[2]-0.8660254037844386*bxby[2]*ux[2])*dx1*wvpar+((-0.8660254037844386*bxbz[0]*uz[3])-0.8660254037844386*bxby[0]*uy[3]-0.8660254037844386*bxbx[0]*ux[3]-0.8660254037844386*uz[1]*bxbz[2]-0.8660254037844386*uy[1]*bxby[2]-0.8660254037844386*ux[1]*bxbx[2])*dx0*wvpar+((-0.4330127018922193*bybz[3]*uz[3])-0.4330127018922193*byby[3]*uy[3]-0.4330127018922193*bxby[3]*ux[3]-0.4330127018922193*bybz[2]*uz[2]-0.4330127018922193*byby[2]*uy[2]-0.4330127018922193*bxby[2]*ux[2])*dvpar*dx1+(0.8660254037844386*Pyz[3]*rho_inv_bz[3]+0.8660254037844386*Pyy[3]*rho_inv_by[3]+0.8660254037844386*Pxy[3]*rho_inv_bx[3]+0.8660254037844386*Pyz[2]*rho_inv_bz[2]+0.8660254037844386*Pyy[2]*rho_inv_by[2]+0.8660254037844386*Pxy[2]*rho_inv_bx[2])*dx1+((-0.4330127018922193*bxbz[0]*uz[3])-0.4330127018922193*bxby[0]*uy[3]-0.4330127018922193*bxbx[0]*ux[3]-0.4330127018922193*uz[1]*bxbz[2]-0.4330127018922193*uy[1]*bxby[2]-0.4330127018922193*ux[1]*bxbx[2])*dvpar*dx0+(0.8660254037844386*rho_inv_bz[0]*Pxz[3]+0.8660254037844386*rho_inv_by[0]*Pxy[3]+0.8660254037844386*rho_inv_bx[0]*Pxx[3]+0.8660254037844386*Pxz[1]*rho_inv_bz[2]+0.8660254037844386*Pxy[1]*rho_inv_by[2]+0.8660254037844386*Pxx[1]*rho_inv_bx[2])*dx0; 
  alphaSurf[3] = ((-0.8660254037844386*bybz[2]*uz[3])-0.8660254037844386*byby[2]*uy[3]-0.8660254037844386*bxby[2]*ux[3]-0.8660254037844386*uz[2]*bybz[3]-0.8660254037844386*uy[2]*byby[3]-0.8660254037844386*ux[2]*bxby[3])*dx1*wvpar+((-0.8660254037844386*bxbz[1]*uz[3])-0.8660254037844386*bxby[1]*uy[3]-0.8660254037844386*bxbx[1]*ux[3]-0.8660254037844386*uz[1]*bxbz[3]-0.8660254037844386*uy[1]*bxby[3]-0.8660254037844386*ux[1]*bxbx[3])*dx0*wvpar+((-0.4330127018922193*bybz[2]*uz[3])-0.4330127018922193*byby[2]*uy[3]-0.4330127018922193*bxby[2]*ux[3]-0.4330127018922193*uz[2]*bybz[3]-0.4330127018922193*uy[2]*byby[3]-0.4330127018922193*ux[2]*bxby[3])*dvpar*dx1+(0.8660254037844386*Pyz[2]*rho_inv_bz[3]+0.8660254037844386*Pyy[2]*rho_inv_by[3]+0.8660254037844386*Pxy[2]*rho_inv_bx[3]+0.8660254037844386*rho_inv_bz[2]*Pyz[3]+0.8660254037844386*rho_inv_by[2]*Pyy[3]+0.8660254037844386*rho_inv_bx[2]*Pxy[3])*dx1+((-0.4330127018922193*bxbz[1]*uz[3])-0.4330127018922193*bxby[1]*uy[3]-0.4330127018922193*bxbx[1]*ux[3]-0.4330127018922193*uz[1]*bxbz[3]-0.4330127018922193*uy[1]*bxby[3]-0.4330127018922193*ux[1]*bxbx[3])*dvpar*dx0+(0.8660254037844386*Pxz[1]*rho_inv_bz[3]+0.8660254037844386*Pxy[1]*rho_inv_by[3]+0.8660254037844386*Pxx[1]*rho_inv_bx[3]+0.8660254037844386*rho_inv_bz[1]*Pxz[3]+0.8660254037844386*rho_inv_by[1]*Pxy[3]+0.8660254037844386*rho_inv_bx[1]*Pxx[3])*dx0; 

  if (0.5*alphaSurf[3]-0.5*(alphaSurf[2]+alphaSurf[1])+0.5*alphaSurf[0] > 0) { 
    fUpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(fEdge); 
  } 
  if ((-0.5*alphaSurf[3])+0.5*alphaSurf[2]-0.5*alphaSurf[1]+0.5*alphaSurf[0] > 0) { 
    fUpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(fEdge); 
  } 
  if (0.5*(alphaSurf[1]+alphaSurf[0])-0.5*(alphaSurf[3]+alphaSurf[2]) > 0) { 
    fUpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(fEdge); 
  } 
  if (0.5*(alphaSurf[3]+alphaSurf[2]+alphaSurf[1]+alphaSurf[0]) > 0) { 
    fUpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alphaSurf[3]*fUpwind[3]+0.5*alphaSurf[2]*fUpwind[2]+0.5*alphaSurf[1]*fUpwind[1]+0.5*alphaSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.5*alphaSurf[2]*fUpwind[3]+0.5*fUpwind[2]*alphaSurf[3]+0.5*alphaSurf[0]*fUpwind[1]+0.5*fUpwind[0]*alphaSurf[1]; 
  Ghat[2] = 0.5*alphaSurf[1]*fUpwind[3]+0.5*fUpwind[1]*alphaSurf[3]+0.5*alphaSurf[0]*fUpwind[2]+0.5*fUpwind[0]*alphaSurf[2]; 
  Ghat[3] = 0.5*alphaSurf[0]*fUpwind[3]+0.5*fUpwind[0]*alphaSurf[3]+0.5*alphaSurf[1]*fUpwind[2]+0.5*fUpwind[1]*alphaSurf[2]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += -0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += -0.7071067811865475*Ghat[2]*dv1par; 
  out[3] += -1.224744871391589*Ghat[0]*dv1par; 
  out[4] += -0.7071067811865475*Ghat[3]*dv1par; 
  out[5] += -1.224744871391589*Ghat[1]*dv1par; 
  out[6] += -1.224744871391589*Ghat[2]*dv1par; 
  out[7] += -1.224744871391589*Ghat[3]*dv1par; 
  out[8] += -1.58113883008419*Ghat[0]*dv1par; 
  out[9] += -1.58113883008419*Ghat[1]*dv1par; 
  out[10] += -1.58113883008419*Ghat[2]*dv1par; 
  out[11] += -1.58113883008419*Ghat[3]*dv1par; 

  } else { 

  alphaSurf[0] = ((-0.8660254037844386*bybz[1]*uz[3])-0.8660254037844386*byby[1]*uy[3]-0.8660254037844386*bxby[1]*ux[3]-0.8660254037844386*bybz[0]*uz[2]-0.8660254037844386*byby[0]*uy[2]-0.8660254037844386*bxby[0]*ux[2])*dx1*wvpar+((-0.8660254037844386*bxbz[2]*uz[3])-0.8660254037844386*bxby[2]*uy[3]-0.8660254037844386*bxbx[2]*ux[3]-0.8660254037844386*bxbz[0]*uz[1]-0.8660254037844386*bxby[0]*uy[1]-0.8660254037844386*bxbx[0]*ux[1])*dx0*wvpar+(0.4330127018922193*bybz[1]*uz[3]+0.4330127018922193*byby[1]*uy[3]+0.4330127018922193*bxby[1]*ux[3]+0.4330127018922193*bybz[0]*uz[2]+0.4330127018922193*byby[0]*uy[2]+0.4330127018922193*bxby[0]*ux[2])*dvpar*dx1+(0.8660254037844386*rho_inv_bz[1]*Pyz[3]+0.8660254037844386*rho_inv_by[1]*Pyy[3]+0.8660254037844386*rho_inv_bx[1]*Pxy[3]+0.8660254037844386*rho_inv_bz[0]*Pyz[2]+0.8660254037844386*rho_inv_by[0]*Pyy[2]+0.8660254037844386*rho_inv_bx[0]*Pxy[2])*dx1+(0.4330127018922193*bxbz[2]*uz[3]+0.4330127018922193*bxby[2]*uy[3]+0.4330127018922193*bxbx[2]*ux[3]+0.4330127018922193*bxbz[0]*uz[1]+0.4330127018922193*bxby[0]*uy[1]+0.4330127018922193*bxbx[0]*ux[1])*dvpar*dx0+(0.8660254037844386*rho_inv_bz[2]*Pxz[3]+0.8660254037844386*rho_inv_by[2]*Pxy[3]+0.8660254037844386*rho_inv_bx[2]*Pxx[3]+0.8660254037844386*rho_inv_bz[0]*Pxz[1]+0.8660254037844386*rho_inv_by[0]*Pxy[1]+0.8660254037844386*rho_inv_bx[0]*Pxx[1])*dx0; 
  alphaSurf[1] = ((-0.8660254037844386*bybz[0]*uz[3])-0.8660254037844386*byby[0]*uy[3]-0.8660254037844386*bxby[0]*ux[3]-0.8660254037844386*bybz[1]*uz[2]-0.8660254037844386*byby[1]*uy[2]-0.8660254037844386*bxby[1]*ux[2])*dx1*wvpar+((-0.8660254037844386*bxbz[3]*uz[3])-0.8660254037844386*bxby[3]*uy[3]-0.8660254037844386*bxbx[3]*ux[3]-0.8660254037844386*bxbz[1]*uz[1]-0.8660254037844386*bxby[1]*uy[1]-0.8660254037844386*bxbx[1]*ux[1])*dx0*wvpar+(0.4330127018922193*bybz[0]*uz[3]+0.4330127018922193*byby[0]*uy[3]+0.4330127018922193*bxby[0]*ux[3]+0.4330127018922193*bybz[1]*uz[2]+0.4330127018922193*byby[1]*uy[2]+0.4330127018922193*bxby[1]*ux[2])*dvpar*dx1+(0.8660254037844386*rho_inv_bz[0]*Pyz[3]+0.8660254037844386*rho_inv_by[0]*Pyy[3]+0.8660254037844386*rho_inv_bx[0]*Pxy[3]+0.8660254037844386*rho_inv_bz[1]*Pyz[2]+0.8660254037844386*rho_inv_by[1]*Pyy[2]+0.8660254037844386*rho_inv_bx[1]*Pxy[2])*dx1+(0.4330127018922193*bxbz[3]*uz[3]+0.4330127018922193*bxby[3]*uy[3]+0.4330127018922193*bxbx[3]*ux[3]+0.4330127018922193*bxbz[1]*uz[1]+0.4330127018922193*bxby[1]*uy[1]+0.4330127018922193*bxbx[1]*ux[1])*dvpar*dx0+(0.8660254037844386*Pxz[3]*rho_inv_bz[3]+0.8660254037844386*Pxy[3]*rho_inv_by[3]+0.8660254037844386*Pxx[3]*rho_inv_bx[3]+0.8660254037844386*Pxz[1]*rho_inv_bz[1]+0.8660254037844386*Pxy[1]*rho_inv_by[1]+0.8660254037844386*Pxx[1]*rho_inv_bx[1])*dx0; 
  alphaSurf[2] = ((-0.8660254037844386*bybz[3]*uz[3])-0.8660254037844386*byby[3]*uy[3]-0.8660254037844386*bxby[3]*ux[3]-0.8660254037844386*bybz[2]*uz[2]-0.8660254037844386*byby[2]*uy[2]-0.8660254037844386*bxby[2]*ux[2])*dx1*wvpar+((-0.8660254037844386*bxbz[0]*uz[3])-0.8660254037844386*bxby[0]*uy[3]-0.8660254037844386*bxbx[0]*ux[3]-0.8660254037844386*uz[1]*bxbz[2]-0.8660254037844386*uy[1]*bxby[2]-0.8660254037844386*ux[1]*bxbx[2])*dx0*wvpar+(0.4330127018922193*bybz[3]*uz[3]+0.4330127018922193*byby[3]*uy[3]+0.4330127018922193*bxby[3]*ux[3]+0.4330127018922193*bybz[2]*uz[2]+0.4330127018922193*byby[2]*uy[2]+0.4330127018922193*bxby[2]*ux[2])*dvpar*dx1+(0.8660254037844386*Pyz[3]*rho_inv_bz[3]+0.8660254037844386*Pyy[3]*rho_inv_by[3]+0.8660254037844386*Pxy[3]*rho_inv_bx[3]+0.8660254037844386*Pyz[2]*rho_inv_bz[2]+0.8660254037844386*Pyy[2]*rho_inv_by[2]+0.8660254037844386*Pxy[2]*rho_inv_bx[2])*dx1+(0.4330127018922193*bxbz[0]*uz[3]+0.4330127018922193*bxby[0]*uy[3]+0.4330127018922193*bxbx[0]*ux[3]+0.4330127018922193*uz[1]*bxbz[2]+0.4330127018922193*uy[1]*bxby[2]+0.4330127018922193*ux[1]*bxbx[2])*dvpar*dx0+(0.8660254037844386*rho_inv_bz[0]*Pxz[3]+0.8660254037844386*rho_inv_by[0]*Pxy[3]+0.8660254037844386*rho_inv_bx[0]*Pxx[3]+0.8660254037844386*Pxz[1]*rho_inv_bz[2]+0.8660254037844386*Pxy[1]*rho_inv_by[2]+0.8660254037844386*Pxx[1]*rho_inv_bx[2])*dx0; 
  alphaSurf[3] = ((-0.8660254037844386*bybz[2]*uz[3])-0.8660254037844386*byby[2]*uy[3]-0.8660254037844386*bxby[2]*ux[3]-0.8660254037844386*uz[2]*bybz[3]-0.8660254037844386*uy[2]*byby[3]-0.8660254037844386*ux[2]*bxby[3])*dx1*wvpar+((-0.8660254037844386*bxbz[1]*uz[3])-0.8660254037844386*bxby[1]*uy[3]-0.8660254037844386*bxbx[1]*ux[3]-0.8660254037844386*uz[1]*bxbz[3]-0.8660254037844386*uy[1]*bxby[3]-0.8660254037844386*ux[1]*bxbx[3])*dx0*wvpar+(0.4330127018922193*bybz[2]*uz[3]+0.4330127018922193*byby[2]*uy[3]+0.4330127018922193*bxby[2]*ux[3]+0.4330127018922193*uz[2]*bybz[3]+0.4330127018922193*uy[2]*byby[3]+0.4330127018922193*ux[2]*bxby[3])*dvpar*dx1+(0.8660254037844386*Pyz[2]*rho_inv_bz[3]+0.8660254037844386*Pyy[2]*rho_inv_by[3]+0.8660254037844386*Pxy[2]*rho_inv_bx[3]+0.8660254037844386*rho_inv_bz[2]*Pyz[3]+0.8660254037844386*rho_inv_by[2]*Pyy[3]+0.8660254037844386*rho_inv_bx[2]*Pxy[3])*dx1+(0.4330127018922193*bxbz[1]*uz[3]+0.4330127018922193*bxby[1]*uy[3]+0.4330127018922193*bxbx[1]*ux[3]+0.4330127018922193*uz[1]*bxbz[3]+0.4330127018922193*uy[1]*bxby[3]+0.4330127018922193*ux[1]*bxbx[3])*dvpar*dx0+(0.8660254037844386*Pxz[1]*rho_inv_bz[3]+0.8660254037844386*Pxy[1]*rho_inv_by[3]+0.8660254037844386*Pxx[1]*rho_inv_bx[3]+0.8660254037844386*rho_inv_bz[1]*Pxz[3]+0.8660254037844386*rho_inv_by[1]*Pxy[3]+0.8660254037844386*rho_inv_bx[1]*Pxx[3])*dx0; 

  if (0.5*alphaSurf[3]-0.5*(alphaSurf[2]+alphaSurf[1])+0.5*alphaSurf[0] > 0) { 
    fUpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(fSkin); 
  } 
  if ((-0.5*alphaSurf[3])+0.5*alphaSurf[2]-0.5*alphaSurf[1]+0.5*alphaSurf[0] > 0) { 
    fUpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(fSkin); 
  } 
  if (0.5*(alphaSurf[1]+alphaSurf[0])-0.5*(alphaSurf[3]+alphaSurf[2]) > 0) { 
    fUpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(fSkin); 
  } 
  if (0.5*(alphaSurf[3]+alphaSurf[2]+alphaSurf[1]+alphaSurf[0]) > 0) { 
    fUpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alphaSurf[3]*fUpwind[3]+0.5*alphaSurf[2]*fUpwind[2]+0.5*alphaSurf[1]*fUpwind[1]+0.5*alphaSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.5*alphaSurf[2]*fUpwind[3]+0.5*fUpwind[2]*alphaSurf[3]+0.5*alphaSurf[0]*fUpwind[1]+0.5*fUpwind[0]*alphaSurf[1]; 
  Ghat[2] = 0.5*alphaSurf[1]*fUpwind[3]+0.5*fUpwind[1]*alphaSurf[3]+0.5*alphaSurf[0]*fUpwind[2]+0.5*fUpwind[0]*alphaSurf[2]; 
  Ghat[3] = 0.5*alphaSurf[0]*fUpwind[3]+0.5*fUpwind[0]*alphaSurf[3]+0.5*alphaSurf[1]*fUpwind[2]+0.5*fUpwind[1]*alphaSurf[2]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += 0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += 0.7071067811865475*Ghat[2]*dv1par; 
  out[3] += -1.224744871391589*Ghat[0]*dv1par; 
  out[4] += 0.7071067811865475*Ghat[3]*dv1par; 
  out[5] += -1.224744871391589*Ghat[1]*dv1par; 
  out[6] += -1.224744871391589*Ghat[2]*dv1par; 
  out[7] += -1.224744871391589*Ghat[3]*dv1par; 
  out[8] += 1.58113883008419*Ghat[0]*dv1par; 
  out[9] += 1.58113883008419*Ghat[1]*dv1par; 
  out[10] += 1.58113883008419*Ghat[2]*dv1par; 
  out[11] += 1.58113883008419*Ghat[3]*dv1par; 

  } 
} 
