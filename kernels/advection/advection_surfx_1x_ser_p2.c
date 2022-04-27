#include <gkyl_advection_kernels.h> 
#include <gkyl_basis_ser_1x_p2_surfx1_eval_quad.h> 
GKYL_CU_DH void advection_surfx_1x_ser_p2(const double *w, const double *dxv, const double *ul, const double *uc, const double *ur, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // ul/uc/ur:  Advection velocity in left/center/right cells.
  // fl/fc/fr:  Input function in left/center/right cells.
  // out:       Incremented function in center cell.
  const double dx1 = 2.0/dxv[0]; 
  const double *ul_0 = &ul[0]; 
  const double *uc_0 = &uc[0]; 
  const double *ur_0 = &ur[0]; 
  double u_l_r = ser_1x_p2_surfx1_eval_quad_node_0_r(ul_0); 
  double u_c_l = ser_1x_p2_surfx1_eval_quad_node_0_l(uc_0); 
  double u_c_r = ser_1x_p2_surfx1_eval_quad_node_0_r(uc_0); 
  double u_r_l = ser_1x_p2_surfx1_eval_quad_node_0_l(ur_0); 
  double u_max_l = fabs(fmax(u_l_r, u_c_l)); 
  double u_max_r= fabs(fmax(u_c_r, u_r_l)); 
  out[0] += dx1*((0.5590169943749475*fr[2]-0.5590169943749475*fc[2]-0.4330127018922193*(fr[1]+fc[1])+0.25*fr[0]-0.25*fc[0])*u_max_r+(0.5590169943749475*fl[2]-0.5590169943749475*fc[2]+0.4330127018922193*(fl[1]+fc[1])+0.25*fl[0]-0.25*fc[0])*u_max_l+((-0.883883476483184*fr[2])+0.6846531968814573*fr[1]-0.3952847075210473*fr[0])*ur_0[2]+(0.883883476483184*fl[2]+0.6846531968814573*fl[1]+0.3952847075210473*fl[0])*ul_0[2]-1.369306393762915*fc[1]*uc_0[2]+(0.6846531968814573*ur_0[1]-0.3952847075210473*ur_0[0])*fr[2]+(0.6846531968814573*ul_0[1]+0.3952847075210473*ul_0[0])*fl[2]-1.369306393762915*uc_0[1]*fc[2]+(0.3061862178478971*fr[0]-0.5303300858899105*fr[1])*ur_0[1]+(0.5303300858899105*fl[1]+0.3061862178478971*fl[0])*ul_0[1]-0.6123724356957944*fc[0]*uc_0[1]+0.3061862178478971*(ur_0[0]*fr[1]+ul_0[0]*fl[1])-0.6123724356957944*uc_0[0]*fc[1]-0.1767766952966368*fr[0]*ur_0[0]+0.1767766952966368*fl[0]*ul_0[0]); 
  out[1] += dx1*((0.9682458365518543*fr[2]-0.9682458365518543*fc[2]-0.75*(fr[1]+fc[1])+0.4330127018922193*fr[0]-0.4330127018922193*fc[0])*u_max_r+((-0.9682458365518543*fl[2])+0.9682458365518543*fc[2]-0.75*(fl[1]+fc[1])-0.4330127018922193*fl[0]+0.4330127018922193*fc[0])*u_max_l+((-1.530931089239486*fr[2])+1.185854122563142*fr[1]-0.6846531968814573*fr[0])*ur_0[2]+((-1.530931089239486*fl[2])-1.185854122563142*fl[1]-0.6846531968814573*fl[0])*ul_0[2]+((-3.061862178478972*fc[2])-1.369306393762915*fc[0])*uc_0[2]+(1.185854122563142*ur_0[1]-0.6846531968814573*ur_0[0])*fr[2]+((-1.185854122563142*ul_0[1])-0.6846531968814573*ul_0[0])*fl[2]-1.369306393762915*uc_0[0]*fc[2]+(0.5303300858899105*fr[0]-0.9185586535436913*fr[1])*ur_0[1]+((-0.9185586535436913*fl[1])-0.5303300858899105*fl[0])*ul_0[1]-1.837117307087383*fc[1]*uc_0[1]+0.5303300858899105*ur_0[0]*fr[1]-0.5303300858899105*ul_0[0]*fl[1]-0.3061862178478971*(fr[0]*ur_0[0]+fl[0]*ul_0[0])-0.6123724356957944*fc[0]*uc_0[0]); 
  out[2] += dx1*((1.25*fr[2]-1.25*fc[2]-0.9682458365518543*(fr[1]+fc[1])+0.5590169943749475*fr[0]-0.5590169943749475*fc[0])*u_max_r+(1.25*fl[2]-1.25*fc[2]+0.9682458365518543*(fl[1]+fc[1])+0.5590169943749475*fl[0]-0.5590169943749475*fc[0])*u_max_l+((-1.976423537605237*fr[2])+1.530931089239486*fr[1]-0.883883476483184*fr[0])*ur_0[2]+(1.976423537605237*fl[2]+1.530931089239486*fl[1]+0.883883476483184*fl[0])*ul_0[2]-3.061862178478972*fc[1]*uc_0[2]+(1.530931089239486*ur_0[1]-0.883883476483184*ur_0[0])*fr[2]+(1.530931089239486*ul_0[1]+0.883883476483184*ul_0[0])*fl[2]-3.061862178478972*uc_0[1]*fc[2]+(0.6846531968814573*fr[0]-1.185854122563142*fr[1])*ur_0[1]+(1.185854122563142*fl[1]+0.6846531968814573*fl[0])*ul_0[1]-1.369306393762915*fc[0]*uc_0[1]+0.6846531968814573*(ur_0[0]*fr[1]+ul_0[0]*fl[1])-1.369306393762915*uc_0[0]*fc[1]-0.3952847075210473*fr[0]*ur_0[0]+0.3952847075210473*fl[0]*ul_0[0]); 

} 
