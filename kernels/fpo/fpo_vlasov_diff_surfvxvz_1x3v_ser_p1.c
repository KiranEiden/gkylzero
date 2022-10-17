#include <gkyl_fpo_vlasov_kernels.h>

GKYL_CU_DH void
fpo_vlasov_diff_surfvxvz_1x3v_ser_p1(const double* w, const double* dx,
  const double* g[], const double* f[], double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // g: 
  // f: 
  // out: Incremented output

  const double Jvxvz = 4/dx[1]/dx[3];

  const double* gllc = g[1];
  const double* glcl = g[3];
  const double* glcc = g[4];
  const double* glcu = g[5];
  const double* gluc = g[7];
  const double* gcll = g[9];
  const double* gclc = g[10];
  const double* gclu = g[11];
  const double* gccl = g[12];
  const double* gccc = g[13];
  const double* gccu = g[14];
  const double* gcul = g[15];
  const double* gcuc = g[16];
  const double* gcuu = g[17];
  const double* gulc = g[19];
  const double* gucl = g[21];
  const double* gucc = g[22];
  const double* gucu = g[23];
  const double* guuc = g[25];

  const double* fllc = f[1];
  const double* flcl = f[3];
  const double* flcc = f[4];
  const double* flcu = f[5];
  const double* fluc = f[7];
  const double* fcll = f[9];
  const double* fclc = f[10];
  const double* fclu = f[11];
  const double* fccl = f[12];
  const double* fccc = f[13];
  const double* fccu = f[14];
  const double* fcul = f[15];
  const double* fcuc = f[16];
  const double* fcuu = f[17];
  const double* fulc = f[19];
  const double* fucl = f[21];
  const double* fucc = f[22];
  const double* fucu = f[23];
  const double* fuuc = f[25];

  double D_proj1_l[8];
  D_proj1_l[0] = (-1.325825214724776*glcc[9])-1.325825214724776*gccc[9]-1.377837980315537*glcc[4]+1.377837980315537*gccc[4];
  D_proj1_l[1] = (-1.325825214724776*glcc[12])-1.325825214724776*gccc[12]-1.377837980315537*glcc[8]+1.377837980315537*gccc[8];
  D_proj1_l[2] = (-1.325825214724776*glcc[14])-1.325825214724776*gccc[14]-1.377837980315537*glcc[10]+1.377837980315537*gccc[10];
  D_proj1_l[3] = 0.0;
  D_proj1_l[4] = (-1.325825214724776*glcc[15])-1.325825214724776*gccc[15]-1.377837980315537*glcc[13]+1.377837980315537*gccc[13];
  D_proj1_l[5] = 0.0;
  D_proj1_l[6] = 0.0;
  D_proj1_l[7] = 0.0;

  double D_proj1_u[8];
  D_proj1_u[0] = (-1.325825214724776*gucc[9])-1.325825214724776*gccc[9]+1.377837980315537*gucc[4]-1.377837980315537*gccc[4];
  D_proj1_u[1] = (-1.325825214724776*gucc[12])-1.325825214724776*gccc[12]+1.377837980315537*gucc[8]-1.377837980315537*gccc[8];
  D_proj1_u[2] = (-1.325825214724776*gucc[14])-1.325825214724776*gccc[14]+1.377837980315537*gucc[10]-1.377837980315537*gccc[10];
  D_proj1_u[3] = 0.0;
  D_proj1_u[4] = (-1.325825214724776*gucc[15])-1.325825214724776*gccc[15]+1.377837980315537*gucc[13]-1.377837980315537*gccc[13];
  D_proj1_u[5] = 0.0;
  D_proj1_u[6] = 0.0;
  D_proj1_u[7] = 0.0;

  double df_proj1_l[8];
  df_proj1_l[0] = (-0.1178511301977579*flcu[9])-0.1178511301977579*flcl[9]+0.2357022603955158*flcc[9]+0.1178511301977579*fccu[9]+0.1178511301977579*fccl[9]-0.2357022603955158*fccc[9]-0.1020620726159657*flcu[4]-0.1020620726159657*flcl[4]+0.2041241452319315*flcc[4]-0.1020620726159657*fccu[4]-0.1020620726159657*fccl[4]+0.2041241452319315*fccc[4]+0.1020620726159657*flcu[2]-0.1020620726159657*flcl[2]-0.1020620726159657*fccu[2]+0.1020620726159657*fccl[2]+0.0883883476483184*flcu[0]-0.0883883476483184*flcl[0]+0.0883883476483184*fccu[0]-0.0883883476483184*fccl[0];
  df_proj1_l[1] = (-0.1178511301977579*flcu[12])-0.1178511301977579*flcl[12]+0.2357022603955158*flcc[12]+0.1178511301977579*fccu[12]+0.1178511301977579*fccl[12]-0.2357022603955158*fccc[12]-0.1020620726159657*flcu[8]-0.1020620726159657*flcl[8]+0.2041241452319315*flcc[8]-0.1020620726159657*fccu[8]-0.1020620726159657*fccl[8]+0.2041241452319315*fccc[8]+0.1020620726159657*flcu[5]-0.1020620726159657*flcl[5]-0.1020620726159657*fccu[5]+0.1020620726159657*fccl[5]+0.0883883476483184*flcu[1]-0.0883883476483184*flcl[1]+0.0883883476483184*fccu[1]-0.0883883476483184*fccl[1];
  df_proj1_l[2] = (-0.1178511301977579*flcu[14])-0.1178511301977579*flcl[14]+0.2357022603955158*flcc[14]+0.1178511301977579*fccu[14]+0.1178511301977579*fccl[14]-0.2357022603955158*fccc[14]-0.1020620726159657*flcu[10]-0.1020620726159657*flcl[10]+0.2041241452319315*flcc[10]-0.1020620726159657*fccu[10]-0.1020620726159657*fccl[10]+0.2041241452319315*fccc[10]+0.1020620726159657*flcu[7]-0.1020620726159657*flcl[7]-0.1020620726159657*fccu[7]+0.1020620726159657*fccl[7]+0.0883883476483184*flcu[3]-0.0883883476483184*flcl[3]+0.0883883476483184*fccu[3]-0.0883883476483184*fccl[3];
  df_proj1_l[3] = (-0.2041241452319315*flcu[9])+0.2041241452319315*flcl[9]+0.2041241452319315*fccu[9]-0.2041241452319315*fccl[9]-0.1767766952966368*flcu[4]+0.1767766952966368*flcl[4]-0.1767766952966368*fccu[4]+0.1767766952966368*fccl[4]+0.1767766952966368*flcu[2]+0.1767766952966368*flcl[2]-0.3535533905932737*flcc[2]-0.1767766952966368*fccu[2]-0.1767766952966368*fccl[2]+0.3535533905932737*fccc[2]+0.1530931089239486*flcu[0]+0.1530931089239486*flcl[0]-0.3061862178478971*flcc[0]+0.1530931089239486*fccu[0]+0.1530931089239486*fccl[0]-0.3061862178478971*fccc[0];
  df_proj1_l[4] = (-0.1178511301977579*flcu[15])-0.1178511301977579*flcl[15]+0.2357022603955158*flcc[15]+0.1178511301977579*fccu[15]+0.1178511301977579*fccl[15]-0.2357022603955158*fccc[15]-0.1020620726159657*flcu[13]-0.1020620726159657*flcl[13]+0.2041241452319315*flcc[13]-0.1020620726159657*fccu[13]-0.1020620726159657*fccl[13]+0.2041241452319315*fccc[13]+0.1020620726159657*flcu[11]-0.1020620726159657*flcl[11]-0.1020620726159657*fccu[11]+0.1020620726159657*fccl[11]+0.0883883476483184*flcu[6]-0.0883883476483184*flcl[6]+0.0883883476483184*fccu[6]-0.0883883476483184*fccl[6];
  df_proj1_l[5] = (-0.2041241452319315*flcu[12])+0.2041241452319315*flcl[12]+0.2041241452319315*fccu[12]-0.2041241452319315*fccl[12]-0.1767766952966368*flcu[8]+0.1767766952966368*flcl[8]-0.1767766952966368*fccu[8]+0.1767766952966368*fccl[8]+0.1767766952966368*flcu[5]+0.1767766952966368*flcl[5]-0.3535533905932737*flcc[5]-0.1767766952966368*fccu[5]-0.1767766952966368*fccl[5]+0.3535533905932737*fccc[5]+0.1530931089239486*flcu[1]+0.1530931089239486*flcl[1]-0.3061862178478971*flcc[1]+0.1530931089239486*fccu[1]+0.1530931089239486*fccl[1]-0.3061862178478971*fccc[1];
  df_proj1_l[6] = (-0.2041241452319315*flcu[14])+0.2041241452319315*flcl[14]+0.2041241452319315*fccu[14]-0.2041241452319315*fccl[14]-0.1767766952966368*flcu[10]+0.1767766952966368*flcl[10]-0.1767766952966368*fccu[10]+0.1767766952966368*fccl[10]+0.1767766952966368*flcu[7]+0.1767766952966368*flcl[7]-0.3535533905932737*flcc[7]-0.1767766952966368*fccu[7]-0.1767766952966368*fccl[7]+0.3535533905932737*fccc[7]+0.1530931089239486*flcu[3]+0.1530931089239486*flcl[3]-0.3061862178478971*flcc[3]+0.1530931089239486*fccu[3]+0.1530931089239486*fccl[3]-0.3061862178478971*fccc[3];
  df_proj1_l[7] = (-0.2041241452319315*flcu[15])+0.2041241452319315*flcl[15]+0.2041241452319315*fccu[15]-0.2041241452319315*fccl[15]-0.1767766952966368*flcu[13]+0.1767766952966368*flcl[13]-0.1767766952966368*fccu[13]+0.1767766952966368*fccl[13]+0.1767766952966368*flcu[11]+0.1767766952966368*flcl[11]-0.3535533905932737*flcc[11]-0.1767766952966368*fccu[11]-0.1767766952966368*fccl[11]+0.3535533905932737*fccc[11]+0.1530931089239486*flcu[6]+0.1530931089239486*flcl[6]-0.3061862178478971*flcc[6]+0.1530931089239486*fccu[6]+0.1530931089239486*fccl[6]-0.3061862178478971*fccc[6];

  double df_proj1_u[8];
  df_proj1_u[0] = 0.1178511301977579*fucu[9]+0.1178511301977579*fucl[9]-0.2357022603955158*fucc[9]-0.1178511301977579*fccu[9]-0.1178511301977579*fccl[9]+0.2357022603955158*fccc[9]-0.1020620726159657*fucu[4]-0.1020620726159657*fucl[4]+0.2041241452319315*fucc[4]-0.1020620726159657*fccu[4]-0.1020620726159657*fccl[4]+0.2041241452319315*fccc[4]-0.1020620726159657*fucu[2]+0.1020620726159657*fucl[2]+0.1020620726159657*fccu[2]-0.1020620726159657*fccl[2]+0.0883883476483184*fucu[0]-0.0883883476483184*fucl[0]+0.0883883476483184*fccu[0]-0.0883883476483184*fccl[0];
  df_proj1_u[1] = 0.1178511301977579*fucu[12]+0.1178511301977579*fucl[12]-0.2357022603955158*fucc[12]-0.1178511301977579*fccu[12]-0.1178511301977579*fccl[12]+0.2357022603955158*fccc[12]-0.1020620726159657*fucu[8]-0.1020620726159657*fucl[8]+0.2041241452319315*fucc[8]-0.1020620726159657*fccu[8]-0.1020620726159657*fccl[8]+0.2041241452319315*fccc[8]-0.1020620726159657*fucu[5]+0.1020620726159657*fucl[5]+0.1020620726159657*fccu[5]-0.1020620726159657*fccl[5]+0.0883883476483184*fucu[1]-0.0883883476483184*fucl[1]+0.0883883476483184*fccu[1]-0.0883883476483184*fccl[1];
  df_proj1_u[2] = 0.1178511301977579*fucu[14]+0.1178511301977579*fucl[14]-0.2357022603955158*fucc[14]-0.1178511301977579*fccu[14]-0.1178511301977579*fccl[14]+0.2357022603955158*fccc[14]-0.1020620726159657*fucu[10]-0.1020620726159657*fucl[10]+0.2041241452319315*fucc[10]-0.1020620726159657*fccu[10]-0.1020620726159657*fccl[10]+0.2041241452319315*fccc[10]-0.1020620726159657*fucu[7]+0.1020620726159657*fucl[7]+0.1020620726159657*fccu[7]-0.1020620726159657*fccl[7]+0.0883883476483184*fucu[3]-0.0883883476483184*fucl[3]+0.0883883476483184*fccu[3]-0.0883883476483184*fccl[3];
  df_proj1_u[3] = 0.2041241452319315*fucu[9]-0.2041241452319315*fucl[9]-0.2041241452319315*fccu[9]+0.2041241452319315*fccl[9]-0.1767766952966368*fucu[4]+0.1767766952966368*fucl[4]-0.1767766952966368*fccu[4]+0.1767766952966368*fccl[4]-0.1767766952966368*fucu[2]-0.1767766952966368*fucl[2]+0.3535533905932737*fucc[2]+0.1767766952966368*fccu[2]+0.1767766952966368*fccl[2]-0.3535533905932737*fccc[2]+0.1530931089239486*fucu[0]+0.1530931089239486*fucl[0]-0.3061862178478971*fucc[0]+0.1530931089239486*fccu[0]+0.1530931089239486*fccl[0]-0.3061862178478971*fccc[0];
  df_proj1_u[4] = 0.1178511301977579*fucu[15]+0.1178511301977579*fucl[15]-0.2357022603955158*fucc[15]-0.1178511301977579*fccu[15]-0.1178511301977579*fccl[15]+0.2357022603955158*fccc[15]-0.1020620726159657*fucu[13]-0.1020620726159657*fucl[13]+0.2041241452319315*fucc[13]-0.1020620726159657*fccu[13]-0.1020620726159657*fccl[13]+0.2041241452319315*fccc[13]-0.1020620726159657*fucu[11]+0.1020620726159657*fucl[11]+0.1020620726159657*fccu[11]-0.1020620726159657*fccl[11]+0.0883883476483184*fucu[6]-0.0883883476483184*fucl[6]+0.0883883476483184*fccu[6]-0.0883883476483184*fccl[6];
  df_proj1_u[5] = 0.2041241452319315*fucu[12]-0.2041241452319315*fucl[12]-0.2041241452319315*fccu[12]+0.2041241452319315*fccl[12]-0.1767766952966368*fucu[8]+0.1767766952966368*fucl[8]-0.1767766952966368*fccu[8]+0.1767766952966368*fccl[8]-0.1767766952966368*fucu[5]-0.1767766952966368*fucl[5]+0.3535533905932737*fucc[5]+0.1767766952966368*fccu[5]+0.1767766952966368*fccl[5]-0.3535533905932737*fccc[5]+0.1530931089239486*fucu[1]+0.1530931089239486*fucl[1]-0.3061862178478971*fucc[1]+0.1530931089239486*fccu[1]+0.1530931089239486*fccl[1]-0.3061862178478971*fccc[1];
  df_proj1_u[6] = 0.2041241452319315*fucu[14]-0.2041241452319315*fucl[14]-0.2041241452319315*fccu[14]+0.2041241452319315*fccl[14]-0.1767766952966368*fucu[10]+0.1767766952966368*fucl[10]-0.1767766952966368*fccu[10]+0.1767766952966368*fccl[10]-0.1767766952966368*fucu[7]-0.1767766952966368*fucl[7]+0.3535533905932737*fucc[7]+0.1767766952966368*fccu[7]+0.1767766952966368*fccl[7]-0.3535533905932737*fccc[7]+0.1530931089239486*fucu[3]+0.1530931089239486*fucl[3]-0.3061862178478971*fucc[3]+0.1530931089239486*fccu[3]+0.1530931089239486*fccl[3]-0.3061862178478971*fccc[3];
  df_proj1_u[7] = 0.2041241452319315*fucu[15]-0.2041241452319315*fucl[15]-0.2041241452319315*fccu[15]+0.2041241452319315*fccl[15]-0.1767766952966368*fucu[13]+0.1767766952966368*fucl[13]-0.1767766952966368*fccu[13]+0.1767766952966368*fccl[13]-0.1767766952966368*fucu[11]-0.1767766952966368*fucl[11]+0.3535533905932737*fucc[11]+0.1767766952966368*fccu[11]+0.1767766952966368*fccl[11]-0.3535533905932737*fccc[11]+0.1530931089239486*fucu[6]+0.1530931089239486*fucl[6]-0.3061862178478971*fucc[6]+0.1530931089239486*fccu[6]+0.1530931089239486*fccl[6]-0.3061862178478971*fccc[6];

  double D_proj2_l[8];
  D_proj2_l[0] = (-1.325825214724776*gccl[9])-1.325825214724776*gccc[9]-1.377837980315537*gccl[2]+1.377837980315537*gccc[2];
  D_proj2_l[1] = (-1.325825214724776*gccl[12])-1.325825214724776*gccc[12]-1.377837980315537*gccl[5]+1.377837980315537*gccc[5];
  D_proj2_l[2] = 0.0;
  D_proj2_l[3] = (-1.325825214724776*gccl[14])-1.325825214724776*gccc[14]-1.377837980315537*gccl[7]+1.377837980315537*gccc[7];
  D_proj2_l[4] = 0.0;
  D_proj2_l[5] = (-1.325825214724776*gccl[15])-1.325825214724776*gccc[15]-1.377837980315537*gccl[11]+1.377837980315537*gccc[11];
  D_proj2_l[6] = 0.0;
  D_proj2_l[7] = 0.0;

  double D_proj2_u[8];
  D_proj2_u[0] = (-1.325825214724776*gccu[9])-1.325825214724776*gccc[9]+1.377837980315537*gccu[2]-1.377837980315537*gccc[2];
  D_proj2_u[1] = (-1.325825214724776*gccu[12])-1.325825214724776*gccc[12]+1.377837980315537*gccu[5]-1.377837980315537*gccc[5];
  D_proj2_u[2] = 0.0;
  D_proj2_u[3] = (-1.325825214724776*gccu[14])-1.325825214724776*gccc[14]+1.377837980315537*gccu[7]-1.377837980315537*gccc[7];
  D_proj2_u[4] = 0.0;
  D_proj2_u[5] = (-1.325825214724776*gccu[15])-1.325825214724776*gccc[15]+1.377837980315537*gccu[11]-1.377837980315537*gccc[11];
  D_proj2_u[6] = 0.0;
  D_proj2_u[7] = 0.0;

  double f_proj2_l[8];
  f_proj2_l[0] = 0.408248290463863*fccl[4]-0.408248290463863*fccc[4]+0.3535533905932737*fccl[0]+0.3535533905932737*fccc[0];
  f_proj2_l[1] = 0.408248290463863*fccl[8]-0.408248290463863*fccc[8]+0.3535533905932737*fccl[1]+0.3535533905932737*fccc[1];
  f_proj2_l[2] = 0.408248290463863*fccl[9]-0.408248290463863*fccc[9]+0.3535533905932737*fccl[2]+0.3535533905932737*fccc[2];
  f_proj2_l[3] = 0.408248290463863*fccl[10]-0.408248290463863*fccc[10]+0.3535533905932737*fccl[3]+0.3535533905932737*fccc[3];
  f_proj2_l[4] = 0.408248290463863*fccl[12]-0.408248290463863*fccc[12]+0.3535533905932737*fccl[5]+0.3535533905932737*fccc[5];
  f_proj2_l[5] = 0.408248290463863*fccl[13]-0.408248290463863*fccc[13]+0.3535533905932737*fccl[6]+0.3535533905932737*fccc[6];
  f_proj2_l[6] = 0.408248290463863*fccl[14]-0.408248290463863*fccc[14]+0.3535533905932737*fccl[7]+0.3535533905932737*fccc[7];
  f_proj2_l[7] = 0.408248290463863*fccl[15]-0.408248290463863*fccc[15]+0.3535533905932737*fccl[11]+0.3535533905932737*fccc[11];

  double f_proj2_u[8];
  f_proj2_u[0] = (-0.408248290463863*fccu[4])+0.408248290463863*fccc[4]+0.3535533905932737*fccu[0]+0.3535533905932737*fccc[0];
  f_proj2_u[1] = (-0.408248290463863*fccu[8])+0.408248290463863*fccc[8]+0.3535533905932737*fccu[1]+0.3535533905932737*fccc[1];
  f_proj2_u[2] = (-0.408248290463863*fccu[9])+0.408248290463863*fccc[9]+0.3535533905932737*fccu[2]+0.3535533905932737*fccc[2];
  f_proj2_u[3] = (-0.408248290463863*fccu[10])+0.408248290463863*fccc[10]+0.3535533905932737*fccu[3]+0.3535533905932737*fccc[3];
  f_proj2_u[4] = (-0.408248290463863*fccu[12])+0.408248290463863*fccc[12]+0.3535533905932737*fccu[5]+0.3535533905932737*fccc[5];
  f_proj2_u[5] = (-0.408248290463863*fccu[13])+0.408248290463863*fccc[13]+0.3535533905932737*fccu[6]+0.3535533905932737*fccc[6];
  f_proj2_u[6] = (-0.408248290463863*fccu[14])+0.408248290463863*fccc[14]+0.3535533905932737*fccu[7]+0.3535533905932737*fccc[7];
  f_proj2_u[7] = (-0.408248290463863*fccu[15])+0.408248290463863*fccc[15]+0.3535533905932737*fccu[11]+0.3535533905932737*fccc[11];

  out[0] +=  Jvxvz*(0.125*D_proj1_u[7]*df_proj1_u[7]-0.125*D_proj1_l[7]*df_proj1_l[7]+0.125*D_proj1_u[6]*df_proj1_u[6]-0.125*D_proj1_l[6]*df_proj1_l[6]+0.125*D_proj1_u[5]*df_proj1_u[5]-0.125*D_proj1_l[5]*df_proj1_l[5]+0.125*D_proj1_u[4]*df_proj1_u[4]-0.125*D_proj1_l[4]*df_proj1_l[4]+0.125*D_proj1_u[3]*df_proj1_u[3]-0.125*D_proj1_l[3]*df_proj1_l[3]+0.125*D_proj1_u[2]*df_proj1_u[2]-0.125*D_proj1_l[2]*df_proj1_l[2]+0.125*D_proj1_u[1]*df_proj1_u[1]-0.125*D_proj1_l[1]*df_proj1_l[1]+0.125*D_proj1_u[0]*df_proj1_u[0]-0.125*D_proj1_l[0]*df_proj1_l[0]);
  out[1] +=  Jvxvz*(0.125*D_proj1_u[6]*df_proj1_u[7]-0.125*D_proj1_l[6]*df_proj1_l[7]+0.125*df_proj1_u[6]*D_proj1_u[7]-0.125*df_proj1_l[6]*D_proj1_l[7]+0.125*D_proj1_u[3]*df_proj1_u[5]-0.125*D_proj1_l[3]*df_proj1_l[5]+0.125*df_proj1_u[3]*D_proj1_u[5]-0.125*df_proj1_l[3]*D_proj1_l[5]+0.125*D_proj1_u[2]*df_proj1_u[4]-0.125*D_proj1_l[2]*df_proj1_l[4]+0.125*df_proj1_u[2]*D_proj1_u[4]-0.125*df_proj1_l[2]*D_proj1_l[4]+0.125*D_proj1_u[0]*df_proj1_u[1]-0.125*D_proj1_l[0]*df_proj1_l[1]+0.125*df_proj1_u[0]*D_proj1_u[1]-0.125*df_proj1_l[0]*D_proj1_l[1]);
  out[2] +=  Jvxvz*((-0.2165063509461096*D_proj2_u[7]*f_proj2_u[7])+0.2165063509461096*D_proj2_l[7]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[7]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[7]*df_proj1_l[7]-0.2165063509461096*D_proj2_u[6]*f_proj2_u[6]+0.2165063509461096*D_proj2_l[6]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[6]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[6]*df_proj1_l[6]-0.2165063509461096*D_proj2_u[5]*f_proj2_u[5]+0.2165063509461096*D_proj2_l[5]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[5]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[5]*df_proj1_l[5]-0.2165063509461096*D_proj2_u[4]*f_proj2_u[4]+0.2165063509461096*D_proj2_l[4]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[4]-0.2165063509461096*D_proj2_u[3]*f_proj2_u[3]+0.2165063509461096*D_proj2_l[3]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[3]-0.2165063509461096*D_proj2_u[2]*f_proj2_u[2]+0.2165063509461096*D_proj2_l[2]*f_proj2_l[2]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[2]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[2]-0.2165063509461096*D_proj2_u[1]*f_proj2_u[1]+0.2165063509461096*D_proj2_l[1]*f_proj2_l[1]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[1]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[1]-0.2165063509461096*D_proj2_u[0]*f_proj2_u[0]+0.2165063509461096*D_proj2_l[0]*f_proj2_l[0]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[0]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[0]);
  out[3] +=  Jvxvz*(0.125*D_proj1_u[5]*df_proj1_u[7]-0.125*D_proj1_l[5]*df_proj1_l[7]+0.125*df_proj1_u[5]*D_proj1_u[7]-0.125*df_proj1_l[5]*D_proj1_l[7]+0.125*D_proj1_u[3]*df_proj1_u[6]-0.125*D_proj1_l[3]*df_proj1_l[6]+0.125*df_proj1_u[3]*D_proj1_u[6]-0.125*df_proj1_l[3]*D_proj1_l[6]+0.125*D_proj1_u[1]*df_proj1_u[4]-0.125*D_proj1_l[1]*df_proj1_l[4]+0.125*df_proj1_u[1]*D_proj1_u[4]-0.125*df_proj1_l[1]*D_proj1_l[4]+0.125*D_proj1_u[0]*df_proj1_u[2]-0.125*D_proj1_l[0]*df_proj1_l[2]+0.125*df_proj1_u[0]*D_proj1_u[2]-0.125*df_proj1_l[0]*D_proj1_l[2]);
  out[4] +=  Jvxvz*(0.125*D_proj1_u[4]*df_proj1_u[7]-0.125*D_proj1_l[4]*df_proj1_l[7]+0.125*df_proj1_u[4]*D_proj1_u[7]-0.125*df_proj1_l[4]*D_proj1_l[7]+0.125*D_proj1_u[2]*df_proj1_u[6]-0.125*D_proj1_l[2]*df_proj1_l[6]+0.125*df_proj1_u[2]*D_proj1_u[6]-0.125*df_proj1_l[2]*D_proj1_l[6]+0.125*D_proj1_u[1]*df_proj1_u[5]-0.125*D_proj1_l[1]*df_proj1_l[5]+0.125*df_proj1_u[1]*D_proj1_u[5]-0.125*df_proj1_l[1]*D_proj1_l[5]+0.125*D_proj1_u[0]*df_proj1_u[3]-0.125*D_proj1_l[0]*df_proj1_l[3]+0.125*df_proj1_u[0]*D_proj1_u[3]-0.125*df_proj1_l[0]*D_proj1_l[3]);
  out[5] +=  Jvxvz*((-0.2165063509461096*D_proj2_u[6]*f_proj2_u[7])+0.2165063509461096*D_proj2_l[6]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[6]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[6]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[6]*D_proj2_u[7]+0.2165063509461096*f_proj2_l[6]*D_proj2_l[7]+0.2165063509461096*df_proj1_u[6]*D_proj1_u[7]+0.2165063509461096*df_proj1_l[6]*D_proj1_l[7]-0.2165063509461096*D_proj2_u[3]*f_proj2_u[5]+0.2165063509461096*D_proj2_l[3]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[3]*D_proj2_u[5]+0.2165063509461096*f_proj2_l[3]*D_proj2_l[5]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[5]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[5]-0.2165063509461096*D_proj2_u[2]*f_proj2_u[4]+0.2165063509461096*D_proj2_l[2]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[4]-0.2165063509461096*f_proj2_u[2]*D_proj2_u[4]+0.2165063509461096*f_proj2_l[2]*D_proj2_l[4]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[4]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[4]-0.2165063509461096*D_proj2_u[0]*f_proj2_u[1]+0.2165063509461096*D_proj2_l[0]*f_proj2_l[1]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[1]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[1]-0.2165063509461096*f_proj2_u[0]*D_proj2_u[1]+0.2165063509461096*f_proj2_l[0]*D_proj2_l[1]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[1]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[1]);
  out[6] +=  Jvxvz*(0.125*D_proj1_u[3]*df_proj1_u[7]-0.125*D_proj1_l[3]*df_proj1_l[7]+0.125*df_proj1_u[3]*D_proj1_u[7]-0.125*df_proj1_l[3]*D_proj1_l[7]+0.125*D_proj1_u[5]*df_proj1_u[6]-0.125*D_proj1_l[5]*df_proj1_l[6]+0.125*df_proj1_u[5]*D_proj1_u[6]-0.125*df_proj1_l[5]*D_proj1_l[6]+0.125*D_proj1_u[0]*df_proj1_u[4]-0.125*D_proj1_l[0]*df_proj1_l[4]+0.125*df_proj1_u[0]*D_proj1_u[4]-0.125*df_proj1_l[0]*D_proj1_l[4]+0.125*D_proj1_u[1]*df_proj1_u[2]-0.125*D_proj1_l[1]*df_proj1_l[2]+0.125*df_proj1_u[1]*D_proj1_u[2]-0.125*df_proj1_l[1]*D_proj1_l[2]);
  out[7] +=  Jvxvz*((-0.2165063509461096*D_proj2_u[4]*f_proj2_u[7])+0.2165063509461096*D_proj2_l[4]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[5]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[5]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[4]*D_proj2_u[7]+0.2165063509461096*f_proj2_l[4]*D_proj2_l[7]+0.2165063509461096*df_proj1_u[5]*D_proj1_u[7]+0.2165063509461096*df_proj1_l[5]*D_proj1_l[7]-0.2165063509461096*D_proj2_u[2]*f_proj2_u[6]+0.2165063509461096*D_proj2_l[2]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[2]*D_proj2_u[6]+0.2165063509461096*f_proj2_l[2]*D_proj2_l[6]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[6]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[6]-0.2165063509461096*D_proj2_u[1]*f_proj2_u[5]+0.2165063509461096*D_proj2_l[1]*f_proj2_l[5]-0.2165063509461096*f_proj2_u[1]*D_proj2_u[5]+0.2165063509461096*f_proj2_l[1]*D_proj2_l[5]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[4]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[4]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[4]-0.2165063509461096*D_proj2_u[0]*f_proj2_u[3]+0.2165063509461096*D_proj2_l[0]*f_proj2_l[3]-0.2165063509461096*f_proj2_u[0]*D_proj2_u[3]+0.2165063509461096*f_proj2_l[0]*D_proj2_l[3]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[2]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[2]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[2]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[2]);
  out[8] +=  Jvxvz*(0.125*D_proj1_u[2]*df_proj1_u[7]-0.125*D_proj1_l[2]*df_proj1_l[7]+0.125*df_proj1_u[2]*D_proj1_u[7]-0.125*df_proj1_l[2]*D_proj1_l[7]+0.125*D_proj1_u[4]*df_proj1_u[6]-0.125*D_proj1_l[4]*df_proj1_l[6]+0.125*df_proj1_u[4]*D_proj1_u[6]-0.125*df_proj1_l[4]*D_proj1_l[6]+0.125*D_proj1_u[0]*df_proj1_u[5]-0.125*D_proj1_l[0]*df_proj1_l[5]+0.125*df_proj1_u[0]*D_proj1_u[5]-0.125*df_proj1_l[0]*D_proj1_l[5]+0.125*D_proj1_u[1]*df_proj1_u[3]-0.125*D_proj1_l[1]*df_proj1_l[3]+0.125*df_proj1_u[1]*D_proj1_u[3]-0.125*df_proj1_l[1]*D_proj1_l[3]);
  out[9] +=  Jvxvz*((-0.375*D_proj2_u[7]*f_proj2_u[7])-0.375*D_proj2_l[7]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[7]+0.2165063509461096*df_proj1_u[4]*D_proj1_u[7]+0.2165063509461096*df_proj1_l[4]*D_proj1_l[7]-0.375*D_proj2_u[6]*f_proj2_u[6]-0.375*D_proj2_l[6]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[6]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[6]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[6]-0.375*D_proj2_u[5]*f_proj2_u[5]-0.375*D_proj2_l[5]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[5]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[5]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[5]-0.375*D_proj2_u[4]*f_proj2_u[4]-0.375*D_proj2_l[4]*f_proj2_l[4]-0.375*D_proj2_u[3]*f_proj2_u[3]-0.375*D_proj2_l[3]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[3]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[3]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[3]-0.375*D_proj2_u[2]*f_proj2_u[2]-0.375*D_proj2_l[2]*f_proj2_l[2]-0.375*D_proj2_u[1]*f_proj2_u[1]-0.375*D_proj2_l[1]*f_proj2_l[1]-0.375*D_proj2_u[0]*f_proj2_u[0]-0.375*D_proj2_l[0]*f_proj2_l[0]);
  out[10] +=  Jvxvz*(0.125*D_proj1_u[1]*df_proj1_u[7]-0.125*D_proj1_l[1]*df_proj1_l[7]+0.125*df_proj1_u[1]*D_proj1_u[7]-0.125*df_proj1_l[1]*D_proj1_l[7]+0.125*D_proj1_u[0]*df_proj1_u[6]-0.125*D_proj1_l[0]*df_proj1_l[6]+0.125*df_proj1_u[0]*D_proj1_u[6]-0.125*df_proj1_l[0]*D_proj1_l[6]+0.125*D_proj1_u[4]*df_proj1_u[5]-0.125*D_proj1_l[4]*df_proj1_l[5]+0.125*df_proj1_u[4]*D_proj1_u[5]-0.125*df_proj1_l[4]*D_proj1_l[5]+0.125*D_proj1_u[2]*df_proj1_u[3]-0.125*D_proj1_l[2]*df_proj1_l[3]+0.125*df_proj1_u[2]*D_proj1_u[3]-0.125*df_proj1_l[2]*D_proj1_l[3]);
  out[11] +=  Jvxvz*((-0.2165063509461096*D_proj2_u[2]*f_proj2_u[7])+0.2165063509461096*D_proj2_l[2]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[2]*D_proj2_u[7]+0.2165063509461096*f_proj2_l[2]*D_proj2_l[7]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[7]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[7]-0.2165063509461096*D_proj2_u[4]*f_proj2_u[6]+0.2165063509461096*D_proj2_l[4]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[5]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[5]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[4]*D_proj2_u[6]+0.2165063509461096*f_proj2_l[4]*D_proj2_l[6]+0.2165063509461096*df_proj1_u[5]*D_proj1_u[6]+0.2165063509461096*df_proj1_l[5]*D_proj1_l[6]-0.2165063509461096*D_proj2_u[0]*f_proj2_u[5]+0.2165063509461096*D_proj2_l[0]*f_proj2_l[5]-0.2165063509461096*f_proj2_u[0]*D_proj2_u[5]+0.2165063509461096*f_proj2_l[0]*D_proj2_l[5]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[4]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[4]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[4]-0.2165063509461096*D_proj2_u[1]*f_proj2_u[3]+0.2165063509461096*D_proj2_l[1]*f_proj2_l[3]-0.2165063509461096*f_proj2_u[1]*D_proj2_u[3]+0.2165063509461096*f_proj2_l[1]*D_proj2_l[3]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[2]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[2]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[2]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[2]);
  out[12] +=  Jvxvz*((-0.375*D_proj2_u[6]*f_proj2_u[7])-0.375*D_proj2_l[6]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[7]-0.375*f_proj2_u[6]*D_proj2_u[7]-0.375*f_proj2_l[6]*D_proj2_l[7]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[7]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[7]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[6]+0.2165063509461096*df_proj1_u[4]*D_proj1_u[6]+0.2165063509461096*df_proj1_l[4]*D_proj1_l[6]-0.375*D_proj2_u[3]*f_proj2_u[5]-0.375*D_proj2_l[3]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[5]-0.375*f_proj2_u[3]*D_proj2_u[5]-0.375*f_proj2_l[3]*D_proj2_l[5]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[5]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[5]-0.375*D_proj2_u[2]*f_proj2_u[4]-0.375*D_proj2_l[2]*f_proj2_l[4]-0.375*f_proj2_u[2]*D_proj2_u[4]-0.375*f_proj2_l[2]*D_proj2_l[4]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[3]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[3]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[3]-0.375*D_proj2_u[0]*f_proj2_u[1]-0.375*D_proj2_l[0]*f_proj2_l[1]-0.375*f_proj2_u[0]*D_proj2_u[1]-0.375*f_proj2_l[0]*D_proj2_l[1]);
  out[13] +=  Jvxvz*(0.125*D_proj1_u[0]*df_proj1_u[7]-0.125*D_proj1_l[0]*df_proj1_l[7]+0.125*df_proj1_u[0]*D_proj1_u[7]-0.125*df_proj1_l[0]*D_proj1_l[7]+0.125*D_proj1_u[1]*df_proj1_u[6]-0.125*D_proj1_l[1]*df_proj1_l[6]+0.125*df_proj1_u[1]*D_proj1_u[6]-0.125*df_proj1_l[1]*D_proj1_l[6]+0.125*D_proj1_u[2]*df_proj1_u[5]-0.125*D_proj1_l[2]*df_proj1_l[5]+0.125*df_proj1_u[2]*D_proj1_u[5]-0.125*df_proj1_l[2]*D_proj1_l[5]+0.125*D_proj1_u[3]*df_proj1_u[4]-0.125*D_proj1_l[3]*df_proj1_l[4]+0.125*df_proj1_u[3]*D_proj1_u[4]-0.125*df_proj1_l[3]*D_proj1_l[4]);
  out[14] +=  Jvxvz*((-0.375*D_proj2_u[4]*f_proj2_u[7])-0.375*D_proj2_l[4]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[7]-0.375*f_proj2_u[4]*D_proj2_u[7]-0.375*f_proj2_l[4]*D_proj2_l[7]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[7]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[7]-0.375*D_proj2_u[2]*f_proj2_u[6]-0.375*D_proj2_l[2]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[6]-0.375*f_proj2_u[2]*D_proj2_u[6]-0.375*f_proj2_l[2]*D_proj2_l[6]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[6]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[6]-0.375*D_proj2_u[1]*f_proj2_u[5]-0.375*D_proj2_l[1]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[5]-0.375*f_proj2_u[1]*D_proj2_u[5]-0.375*f_proj2_l[1]*D_proj2_l[5]+0.2165063509461096*df_proj1_u[4]*D_proj1_u[5]+0.2165063509461096*df_proj1_l[4]*D_proj1_l[5]-0.375*D_proj2_u[0]*f_proj2_u[3]-0.375*D_proj2_l[0]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[3]-0.375*f_proj2_u[0]*D_proj2_u[3]-0.375*f_proj2_l[0]*D_proj2_l[3]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[3]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[3]);
  out[15] +=  Jvxvz*((-0.375*D_proj2_u[2]*f_proj2_u[7])-0.375*D_proj2_l[2]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[7]-0.375*f_proj2_u[2]*D_proj2_u[7]-0.375*f_proj2_l[2]*D_proj2_l[7]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[7]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[7]-0.375*D_proj2_u[4]*f_proj2_u[6]-0.375*D_proj2_l[4]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[6]-0.375*f_proj2_u[4]*D_proj2_u[6]-0.375*f_proj2_l[4]*D_proj2_l[6]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[6]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[6]-0.375*D_proj2_u[0]*f_proj2_u[5]-0.375*D_proj2_l[0]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[5]-0.375*f_proj2_u[0]*D_proj2_u[5]-0.375*f_proj2_l[0]*D_proj2_l[5]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[5]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[5]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[4]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[4]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[4]-0.375*D_proj2_u[1]*f_proj2_u[3]-0.375*D_proj2_l[1]*f_proj2_l[3]-0.375*f_proj2_u[1]*D_proj2_u[3]-0.375*f_proj2_l[1]*D_proj2_l[3]);
}

