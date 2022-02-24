// Thu Feb 24 09:26:57 2022
#include <gkyl_basis_gk_hyb_kernels.h>
GKYL_CU_DH
void
eval_2d_gk_hyb_p1(const double *z, double *b )
{
  const double z0 = z[0];
  const double z1 = z[1];
  b[0] = 5.0000000000000000e-01;
  b[1] = 8.6602540378443860e-01*z0;
  b[2] = 8.6602540378443860e-01*z1;
  b[3] = 1.5000000000000000e+00*z0*z1;
  b[4] =  1.6770509831248424e+00*(z1*z1)-5.5901699437494745e-01;
  b[5] =  2.9047375096555625e+00*z0*(z1*z1)+-9.6824583655185426e-01*z0;
}

GKYL_CU_DH
double
eval_expand_2d_gk_hyb_p1(const double *z, const double *f )
{
  const double z0 = z[0];
  const double z1 = z[1];
 return  1.5000000000000000e+00*f[3]*z0*z1+1.6770509831248424e+00*(z1*z1)*f[4]+5.0000000000000000e-01*f[0]+8.6602540378443860e-01*f[2]*z1+-9.6824583655185426e-01*z0*f[5]+2.9047375096555625e+00*z0*(z1*z1)*f[5]+-5.5901699437494745e-01*f[4]+8.6602540378443860e-01*f[1]*z0;
}

GKYL_CU_DH
double
eval_grad_expand_2d_gk_hyb_p1(int dir, const double *z, const double *f )
{
  const double z0 = z[0];
  const double z1 = z[1];
  if (dir == 0)
    return  -9.6824583655185426e-01*f[5]+2.9047375096555625e+00*f[5]*(z1*z1)+1.5000000000000000e+00*z1*f[3]+8.6602540378443860e-01*f[1];

  if (dir == 1)
    return  1.5000000000000000e+00*z0*f[3]+3.3541019662496847e+00*f[4]*z1+8.6602540378443860e-01*f[2]+5.8094750193111251e+00*z0*f[5]*z1;

  return 0.0; // can't happen, suppresses warning

}

GKYL_CU_DH
void
eval_3d_gk_hyb_p1(const double *z, double *b )
{
  const double z0 = z[0];
  const double z1 = z[1];
  const double z2 = z[2];
  b[0] = 3.5355339059327379e-01;
  b[1] = 6.1237243569579447e-01*z0;
  b[2] = 6.1237243569579447e-01*z1;
  b[3] = 6.1237243569579447e-01*z2;
  b[4] = 1.0606601717798212e+00*z0*z1;
  b[5] = 1.0606601717798212e+00*z0*z2;
  b[6] = 1.0606601717798212e+00*z1*z2;
  b[7] = 1.8371173070873836e+00*z0*z1*z2;
  b[8] =  1.1858541225631423e+00*(z1*z1)-3.9528470752104744e-01;
  b[9] =  2.0539595906443728e+00*z0*(z1*z1)+-6.8465319688145765e-01*z0;
  b[10] =  2.0539595906443728e+00*(z1*z1)*z2+-6.8465319688145765e-01*z2;
  b[11] =  -1.1858541225631423e+00*z0*z2+3.5575623676894268e+00*z0*(z1*z1)*z2;
}

GKYL_CU_DH
double
eval_expand_3d_gk_hyb_p1(const double *z, const double *f )
{
  const double z0 = z[0];
  const double z1 = z[1];
  const double z2 = z[2];
 return  1.8371173070873836e+00*f[7]*z0*z1*z2+-6.8465319688145765e-01*f[10]*z2+1.0606601717798212e+00*f[6]*z1*z2+3.5575623676894268e+00*z0*(z1*z1)*z2*f[11]+1.1858541225631423e+00*f[8]*(z1*z1)+6.1237243569579447e-01*z0*f[1]+-6.8465319688145765e-01*z0*f[9]+1.0606601717798212e+00*f[4]*z0*z1+-3.9528470752104744e-01*f[8]+3.5355339059327379e-01*f[0]+-1.1858541225631423e+00*z0*z2*f[11]+2.0539595906443728e+00*f[10]*(z1*z1)*z2+1.0606601717798212e+00*f[5]*z0*z2+6.1237243569579447e-01*f[3]*z2+6.1237243569579447e-01*z1*f[2]+2.0539595906443728e+00*z0*(z1*z1)*f[9];
}

GKYL_CU_DH
double
eval_grad_expand_3d_gk_hyb_p1(int dir, const double *z, const double *f )
{
  const double z0 = z[0];
  const double z1 = z[1];
  const double z2 = z[2];
  if (dir == 0)
    return  1.8371173070873836e+00*z1*f[7]*z2+6.1237243569579447e-01*f[1]+2.0539595906443728e+00*f[9]*(z1*z1)+3.5575623676894268e+00*f[11]*(z1*z1)*z2+-6.8465319688145765e-01*f[9]+1.0606601717798212e+00*f[4]*z1+1.0606601717798212e+00*f[5]*z2+-1.1858541225631423e+00*f[11]*z2;

  if (dir == 1)
    return  1.0606601717798212e+00*f[6]*z2+6.1237243569579447e-01*f[2]+1.0606601717798212e+00*z0*f[4]+7.1151247353788536e+00*z0*f[11]*z1*z2+2.3717082451262845e+00*f[8]*z1+1.8371173070873836e+00*z0*f[7]*z2+4.1079191812887457e+00*z1*f[10]*z2+4.1079191812887457e+00*f[9]*z0*z1;

  if (dir == 2)
    return  1.0606601717798212e+00*z1*f[6]+6.1237243569579447e-01*f[3]+2.0539595906443728e+00*(z1*z1)*f[10]+1.0606601717798212e+00*z0*f[5]+3.5575623676894268e+00*z0*f[11]*(z1*z1)+1.8371173070873836e+00*z0*z1*f[7]+-6.8465319688145765e-01*f[10]+-1.1858541225631423e+00*z0*f[11];

  return 0.0; // can't happen, suppresses warning

}

GKYL_CU_DH
void
eval_4d_gk_hyb_p1(const double *z, double *b )
{
  const double z0 = z[0];
  const double z1 = z[1];
  const double z2 = z[2];
  const double z3 = z[3];
  b[0] = 2.5000000000000000e-01;
  b[1] = 4.3301270189221930e-01*z0;
  b[2] = 4.3301270189221930e-01*z1;
  b[3] = 4.3301270189221930e-01*z2;
  b[4] = 4.3301270189221930e-01*z3;
  b[5] = 7.5000000000000000e-01*z0*z1;
  b[6] = 7.5000000000000000e-01*z0*z2;
  b[7] = 7.5000000000000000e-01*z1*z2;
  b[8] = 7.5000000000000000e-01*z3*z0;
  b[9] = 7.5000000000000000e-01*z3*z1;
  b[10] = 7.5000000000000000e-01*z3*z2;
  b[11] = 1.2990381056766580e+00*z0*z1*z2;
  b[12] = 1.2990381056766580e+00*z3*z0*z1;
  b[13] = 1.2990381056766580e+00*z3*z0*z2;
  b[14] = 1.2990381056766580e+00*z3*z1*z2;
  b[15] = 2.2500000000000000e+00*z3*z0*z1*z2;
  b[16] =  8.3852549156242118e-01*(z2*z2)-2.7950849718747373e-01;
  b[17] =  -4.8412291827592713e-01*z0+1.4523687548277813e+00*z0*(z2*z2);
  b[18] =  -4.8412291827592713e-01*z1+1.4523687548277813e+00*z1*(z2*z2);
  b[19] =  1.4523687548277813e+00*z3*(z2*z2)+-4.8412291827592713e-01*z3;
  b[20] =  2.5155764746872635e+00*z0*z1*(z2*z2)+-8.3852549156242118e-01*z0*z1;
  b[21] =  -8.3852549156242118e-01*z3*z0+2.5155764746872635e+00*z3*z0*(z2*z2);
  b[22] =  -8.3852549156242118e-01*z3*z1+2.5155764746872635e+00*z3*z1*(z2*z2);
  b[23] =  4.3571062644833436e+00*z3*z0*z1*(z2*z2)+-1.4523687548277813e+00*z3*z0*z1;
}

GKYL_CU_DH
double
eval_expand_4d_gk_hyb_p1(const double *z, const double *f )
{
  const double z0 = z[0];
  const double z1 = z[1];
  const double z2 = z[2];
  const double z3 = z[3];
 return  -4.8412291827592713e-01*z3*f[19]+-4.8412291827592713e-01*f[18]*z1+-1.4523687548277813e+00*f[23]*z3*z0*z1+2.2500000000000000e+00*z3*z0*f[15]*z1*z2+2.5155764746872635e+00*z3*z0*f[21]*(z2*z2)+1.2990381056766580e+00*z3*z0*f[13]*z2+1.4523687548277813e+00*z0*f[17]*(z2*z2)+7.5000000000000000e-01*f[6]*z0*z2+2.5155764746872635e+00*z0*z1*f[20]*(z2*z2)+1.2990381056766580e+00*z3*f[12]*z0*z1+4.3301270189221930e-01*z3*f[4]+7.5000000000000000e-01*z0*f[5]*z1+7.5000000000000000e-01*z3*z1*f[9]+2.5155764746872635e+00*z3*f[22]*z1*(z2*z2)+4.3301270189221930e-01*z1*f[2]+1.4523687548277813e+00*f[18]*z1*(z2*z2)+1.2990381056766580e+00*z0*z1*f[11]*z2+-2.7950849718747373e-01*f[16]+1.4523687548277813e+00*z3*f[19]*(z2*z2)+4.3571062644833436e+00*f[23]*z3*z0*z1*(z2*z2)+1.2990381056766580e+00*f[14]*z3*z1*z2+7.5000000000000000e-01*f[10]*z3*z2+-8.3852549156242118e-01*z3*z0*f[21]+7.5000000000000000e-01*z3*z0*f[8]+4.3301270189221930e-01*z0*f[1]+2.5000000000000000e-01*f[0]+8.3852549156242118e-01*f[16]*(z2*z2)+7.5000000000000000e-01*f[7]*z1*z2+-8.3852549156242118e-01*z0*z1*f[20]+-4.8412291827592713e-01*z0*f[17]+-8.3852549156242118e-01*z3*f[22]*z1+4.3301270189221930e-01*f[3]*z2;
}

GKYL_CU_DH
double
eval_grad_expand_4d_gk_hyb_p1(int dir, const double *z, const double *f )
{
  const double z0 = z[0];
  const double z1 = z[1];
  const double z2 = z[2];
  const double z3 = z[3];
  if (dir == 0)
    return  2.5155764746872635e+00*f[20]*z1*(z2*z2)+1.2990381056766580e+00*f[13]*z3*z2+4.3571062644833436e+00*z3*z1*f[23]*(z2*z2)+-8.3852549156242118e-01*z3*f[21]+1.2990381056766580e+00*f[11]*z1*z2+7.5000000000000000e-01*z3*f[8]+-8.3852549156242118e-01*f[20]*z1+4.3301270189221930e-01*f[1]+-1.4523687548277813e+00*z3*z1*f[23]+1.2990381056766580e+00*z3*z1*f[12]+2.2500000000000000e+00*f[15]*z3*z1*z2+7.5000000000000000e-01*z1*f[5]+-4.8412291827592713e-01*f[17]+1.4523687548277813e+00*f[17]*(z2*z2)+2.5155764746872635e+00*z3*f[21]*(z2*z2)+7.5000000000000000e-01*f[6]*z2;

  if (dir == 1)
    return  -1.4523687548277813e+00*z3*z0*f[23]+1.2990381056766580e+00*z3*z0*f[12]+4.3301270189221930e-01*f[2]+2.5155764746872635e+00*f[20]*z0*(z2*z2)+-8.3852549156242118e-01*z3*f[22]+1.2990381056766580e+00*f[11]*z0*z2+7.5000000000000000e-01*f[7]*z2+7.5000000000000000e-01*z0*f[5]+2.2500000000000000e+00*f[15]*z3*z0*z2+-4.8412291827592713e-01*f[18]+4.3571062644833436e+00*z3*z0*f[23]*(z2*z2)+7.5000000000000000e-01*z3*f[9]+2.5155764746872635e+00*z3*(z2*z2)*f[22]+1.2990381056766580e+00*z3*f[14]*z2+-8.3852549156242118e-01*f[20]*z0+1.4523687548277813e+00*f[18]*(z2*z2);

  if (dir == 2)
    return  2.9047375096555625e+00*f[18]*z1*z2+4.3301270189221930e-01*f[3]+5.0311529493745271e+00*z3*z1*z2*f[22]+1.2990381056766580e+00*z3*z1*f[14]+8.7142125289666872e+00*z3*z0*z1*f[23]*z2+5.0311529493745271e+00*z3*z0*f[21]*z2+2.9047375096555625e+00*z0*f[17]*z2+1.2990381056766580e+00*f[13]*z3*z0+5.0311529493745271e+00*f[20]*z0*z1*z2+2.2500000000000000e+00*f[15]*z3*z0*z1+1.2990381056766580e+00*f[11]*z0*z1+7.5000000000000000e-01*z1*f[7]+2.9047375096555625e+00*z3*f[19]*z2+7.5000000000000000e-01*z3*f[10]+7.5000000000000000e-01*z0*f[6]+1.6770509831248424e+00*f[16]*z2;

  if (dir == 3)
    return  2.5155764746872635e+00*z0*f[21]*(z2*z2)+2.2500000000000000e+00*f[15]*z0*z1*z2+-4.8412291827592713e-01*f[19]+7.5000000000000000e-01*f[10]*z2+1.2990381056766580e+00*z1*f[14]*z2+1.2990381056766580e+00*z0*z1*f[12]+-1.4523687548277813e+00*z0*z1*f[23]+7.5000000000000000e-01*f[9]*z1+2.5155764746872635e+00*z1*(z2*z2)*f[22]+-8.3852549156242118e-01*z0*f[21]+1.4523687548277813e+00*f[19]*(z2*z2)+4.3571062644833436e+00*z0*z1*f[23]*(z2*z2)+1.2990381056766580e+00*f[13]*z0*z2+7.5000000000000000e-01*f[8]*z0+-8.3852549156242118e-01*z1*f[22]+4.3301270189221930e-01*f[4];

  return 0.0; // can't happen, suppresses warning

}

GKYL_CU_DH
void
eval_5d_gk_hyb_p1(const double *z, double *b )
{
  const double z0 = z[0];
  const double z1 = z[1];
  const double z2 = z[2];
  const double z3 = z[3];
  const double z4 = z[4];
  b[0] = 1.7677669529663689e-01;
  b[1] = 3.0618621784789724e-01*z0;
  b[2] = 3.0618621784789724e-01*z1;
  b[3] = 3.0618621784789724e-01*z2;
  b[4] = 3.0618621784789724e-01*z3;
  b[5] = 3.0618621784789724e-01*z4;
  b[6] = 5.3033008588991060e-01*z0*z1;
  b[7] = 5.3033008588991060e-01*z0*z2;
  b[8] = 5.3033008588991060e-01*z1*z2;
  b[9] = 5.3033008588991060e-01*z3*z0;
  b[10] = 5.3033008588991060e-01*z3*z1;
  b[11] = 5.3033008588991060e-01*z3*z2;
  b[12] = 5.3033008588991060e-01*z0*z4;
  b[13] = 5.3033008588991060e-01*z4*z1;
  b[14] = 5.3033008588991060e-01*z4*z2;
  b[15] = 5.3033008588991060e-01*z3*z4;
  b[16] = 9.1855865354369182e-01*z0*z1*z2;
  b[17] = 9.1855865354369182e-01*z3*z0*z1;
  b[18] = 9.1855865354369182e-01*z3*z0*z2;
  b[19] = 9.1855865354369182e-01*z3*z1*z2;
  b[20] = 9.1855865354369182e-01*z0*z4*z1;
  b[21] = 9.1855865354369182e-01*z0*z4*z2;
  b[22] = 9.1855865354369182e-01*z4*z1*z2;
  b[23] = 9.1855865354369182e-01*z3*z0*z4;
  b[24] = 9.1855865354369182e-01*z3*z4*z1;
  b[25] = 9.1855865354369182e-01*z3*z4*z2;
  b[26] = 1.5909902576697319e+00*z3*z0*z1*z2;
  b[27] = 1.5909902576697319e+00*z0*z4*z1*z2;
  b[28] = 1.5909902576697319e+00*z3*z0*z4*z1;
  b[29] = 1.5909902576697319e+00*z3*z0*z4*z2;
  b[30] = 1.5909902576697319e+00*z3*z4*z1*z2;
  b[31] = 2.7556759606310752e+00*z3*z0*z4*z1*z2;
  b[32] =  5.9292706128157113e-01*(z3*z3)-1.9764235376052372e-01;
  b[33] =  -3.4232659844072882e-01*z0+1.0269797953221864e+00*(z3*z3)*z0;
  b[34] =  1.0269797953221864e+00*(z3*z3)*z1+-3.4232659844072882e-01*z1;
  b[35] =  1.0269797953221864e+00*(z3*z3)*z2+-3.4232659844072882e-01*z2;
  b[36] =  1.0269797953221864e+00*(z3*z3)*z4+-3.4232659844072882e-01*z4;
  b[37] =  -5.9292706128157113e-01*z0*z1+1.7787811838447134e+00*(z3*z3)*z0*z1;
  b[38] =  1.7787811838447134e+00*(z3*z3)*z0*z2+-5.9292706128157113e-01*z0*z2;
  b[39] =  1.7787811838447134e+00*(z3*z3)*z1*z2+-5.9292706128157113e-01*z1*z2;
  b[40] =  -5.9292706128157113e-01*z0*z4+1.7787811838447134e+00*(z3*z3)*z0*z4;
  b[41] =  -5.9292706128157113e-01*z4*z1+1.7787811838447134e+00*(z3*z3)*z4*z1;
  b[42] =  1.7787811838447134e+00*(z3*z3)*z4*z2+-5.9292706128157113e-01*z4*z2;
  b[43] =  -1.0269797953221864e+00*z0*z1*z2+3.0809393859665595e+00*(z3*z3)*z0*z1*z2;
  b[44] =  3.0809393859665595e+00*(z3*z3)*z0*z4*z1+-1.0269797953221864e+00*z0*z4*z1;
  b[45] =  -1.0269797953221864e+00*z0*z4*z2+3.0809393859665595e+00*(z3*z3)*z0*z4*z2;
  b[46] =  -1.0269797953221864e+00*z4*z1*z2+3.0809393859665595e+00*(z3*z3)*z4*z1*z2;
  b[47] =  5.3363435515341404e+00*(z3*z3)*z0*z4*z1*z2+-1.7787811838447134e+00*z0*z4*z1*z2;
}

GKYL_CU_DH
double
eval_expand_5d_gk_hyb_p1(const double *z, const double *f )
{
  const double z0 = z[0];
  const double z1 = z[1];
  const double z2 = z[2];
  const double z3 = z[3];
  const double z4 = z[4];
 return  9.1855865354369182e-01*f[18]*z3*z0*z2+9.1855865354369182e-01*z3*z4*f[25]*z2+1.5909902576697319e+00*z3*z0*f[26]*z1*z2+3.0618621784789724e-01*f[3]*z2+-5.9292706128157113e-01*z4*z1*f[41]+-1.0269797953221864e+00*z4*z1*f[46]*z2+1.0269797953221864e+00*(z3*z3)*z0*f[33]+1.5909902576697319e+00*z3*z0*z4*z1*f[28]+3.0618621784789724e-01*z3*f[4]+1.5909902576697319e+00*z3*z0*z4*f[29]*z2+9.1855865354369182e-01*z3*z1*f[19]*z2+-1.7787811838447134e+00*f[47]*z0*z4*z1*z2+3.0809393859665595e+00*(z3*z3)*z0*f[43]*z1*z2+1.7787811838447134e+00*(z3*z3)*z0*f[37]*z1+5.3033008588991060e-01*z3*z4*f[15]+1.0269797953221864e+00*(z3*z3)*f[35]*z2+1.0269797953221864e+00*(z3*z3)*f[34]*z1+-1.9764235376052372e-01*f[32]+3.0618621784789724e-01*f[2]*z1+9.1855865354369182e-01*z3*z4*z1*f[24]+-5.9292706128157113e-01*f[42]*z4*z2+5.3033008588991060e-01*z0*f[6]*z1+5.3033008588991060e-01*z3*z0*f[9]+9.1855865354369182e-01*z0*z4*z1*f[20]+3.0809393859665595e+00*(z3*z3)*f[45]*z0*z4*z2+1.0269797953221864e+00*(z3*z3)*z4*f[36]+1.7787811838447134e+00*f[39]*(z3*z3)*z1*z2+9.1855865354369182e-01*z3*z0*f[17]*z1+1.7787811838447134e+00*(z3*z3)*z0*f[38]*z2+5.3033008588991060e-01*z1*f[8]*z2+2.7556759606310752e+00*z3*f[31]*z0*z4*z1*z2+3.0809393859665595e+00*(z3*z3)*z0*z4*z1*f[44]+3.0618621784789724e-01*f[5]*z4+-5.9292706128157113e-01*z0*f[40]*z4+9.1855865354369182e-01*z0*z1*f[16]*z2+5.3363435515341404e+00*f[47]*(z3*z3)*z0*z4*z1*z2+9.1855865354369182e-01*z0*f[21]*z4*z2+-3.4232659844072882e-01*z4*f[36]+-1.0269797953221864e+00*z0*f[43]*z1*z2+1.5909902576697319e+00*f[27]*z0*z4*z1*z2+1.7787811838447134e+00*(z3*z3)*z0*f[40]*z4+5.3033008588991060e-01*f[10]*z3*z1+3.0618621784789724e-01*z0*f[1]+9.1855865354369182e-01*f[22]*z4*z1*z2+-5.9292706128157113e-01*z0*f[38]*z2+-3.4232659844072882e-01*f[34]*z1+3.0809393859665595e+00*(z3*z3)*z4*z1*f[46]*z2+1.7787811838447134e+00*f[42]*(z3*z3)*z4*z2+5.3033008588991060e-01*z3*f[11]*z2+-1.0269797953221864e+00*z0*z4*z1*f[44]+5.9292706128157113e-01*(z3*z3)*f[32]+5.3033008588991060e-01*z0*f[7]*z2+1.7787811838447134e+00*(z3*z3)*z4*z1*f[41]+-1.0269797953221864e+00*f[45]*z0*z4*z2+5.3033008588991060e-01*z4*z1*f[13]+-3.4232659844072882e-01*f[35]*z2+1.7677669529663689e-01*f[0]+-5.9292706128157113e-01*f[39]*z1*z2+-5.9292706128157113e-01*z0*f[37]*z1+9.1855865354369182e-01*f[23]*z3*z0*z4+5.3033008588991060e-01*f[14]*z4*z2+1.5909902576697319e+00*z3*f[30]*z4*z1*z2+-3.4232659844072882e-01*z0*f[33]+5.3033008588991060e-01*f[12]*z0*z4;
}

GKYL_CU_DH
double
eval_grad_expand_5d_gk_hyb_p1(int dir, const double *z, const double *f )
{
  const double z0 = z[0];
  const double z1 = z[1];
  const double z2 = z[2];
  const double z3 = z[3];
  const double z4 = z[4];
  if (dir == 0)
    return  1.5909902576697319e+00*f[26]*z3*z1*z2+-1.7787811838447134e+00*f[47]*z4*z1*z2+5.3033008588991060e-01*f[6]*z1+3.0809393859665595e+00*(z3*z3)*z4*f[44]*z1+1.7787811838447134e+00*(z3*z3)*f[38]*z2+-5.9292706128157113e-01*f[40]*z4+9.1855865354369182e-01*z4*z1*f[20]+2.7556759606310752e+00*z3*f[31]*z4*z1*z2+1.0269797953221864e+00*(z3*z3)*f[33]+9.1855865354369182e-01*f[21]*z4*z2+1.7787811838447134e+00*f[40]*(z3*z3)*z4+5.3033008588991060e-01*z3*f[9]+1.7787811838447134e+00*(z3*z3)*f[37]*z1+-1.0269797953221864e+00*f[45]*z4*z2+9.1855865354369182e-01*f[17]*z3*z1+-1.0269797953221864e+00*z1*z2*f[43]+5.3033008588991060e-01*f[7]*z2+3.0618621784789724e-01*f[1]+9.1855865354369182e-01*z3*f[18]*z2+1.5909902576697319e+00*z3*z4*z1*f[28]+3.0809393859665595e+00*(z3*z3)*z1*z2*f[43]+9.1855865354369182e-01*z3*f[23]*z4+-5.9292706128157113e-01*f[37]*z1+9.1855865354369182e-01*z1*f[16]*z2+5.3363435515341404e+00*(z3*z3)*f[47]*z4*z1*z2+-3.4232659844072882e-01*f[33]+-1.0269797953221864e+00*z4*f[44]*z1+1.5909902576697319e+00*f[27]*z4*z1*z2+1.5909902576697319e+00*z3*z4*f[29]*z2+-5.9292706128157113e-01*f[38]*z2+5.3033008588991060e-01*f[12]*z4+3.0809393859665595e+00*f[45]*(z3*z3)*z4*z2;

  if (dir == 1)
    return  5.3363435515341404e+00*(z3*z3)*z0*f[47]*z4*z2+-5.9292706128157113e-01*f[37]*z0+-5.9292706128157113e-01*z4*f[41]+9.1855865354369182e-01*z0*z4*f[20]+-1.0269797953221864e+00*z4*f[46]*z2+3.0809393859665595e+00*(z3*z3)*z0*z2*f[43]+9.1855865354369182e-01*z3*z4*f[24]+9.1855865354369182e-01*z0*f[16]*z2+1.5909902576697319e+00*f[26]*z3*z0*z2+-3.4232659844072882e-01*f[34]+5.3033008588991060e-01*z3*f[10]+3.0809393859665595e+00*(z3*z3)*z0*z4*f[44]+-5.9292706128157113e-01*f[39]*z2+9.1855865354369182e-01*z3*f[19]*z2+1.7787811838447134e+00*(z3*z3)*f[39]*z2+1.5909902576697319e+00*f[30]*z3*z4*z2+5.3033008588991060e-01*f[8]*z2+1.5909902576697319e+00*z3*z0*z4*f[28]+5.3033008588991060e-01*f[6]*z0+-1.0269797953221864e+00*z0*z4*f[44]+-1.7787811838447134e+00*z0*f[47]*z4*z2+5.3033008588991060e-01*z4*f[13]+1.5909902576697319e+00*z0*f[27]*z4*z2+9.1855865354369182e-01*f[22]*z4*z2+9.1855865354369182e-01*f[17]*z3*z0+2.7556759606310752e+00*z3*z0*f[31]*z4*z2+3.0618621784789724e-01*f[2]+1.0269797953221864e+00*(z3*z3)*f[34]+-1.0269797953221864e+00*z0*z2*f[43]+3.0809393859665595e+00*(z3*z3)*z4*f[46]*z2+1.7787811838447134e+00*(z3*z3)*z4*f[41]+1.7787811838447134e+00*(z3*z3)*f[37]*z0;

  if (dir == 2)
    return  3.0809393859665595e+00*f[45]*(z3*z3)*z0*z4+1.7787811838447134e+00*(z3*z3)*z0*f[38]+-1.0269797953221864e+00*z0*z1*f[43]+5.3033008588991060e-01*f[14]*z4+2.7556759606310752e+00*z3*z0*f[31]*z4*z1+1.5909902576697319e+00*z3*z0*z4*f[29]+9.1855865354369182e-01*f[22]*z4*z1+1.5909902576697319e+00*z0*f[27]*z4*z1+-1.0269797953221864e+00*z4*z1*f[46]+-1.7787811838447134e+00*z0*f[47]*z4*z1+-5.9292706128157113e-01*f[42]*z4+5.3033008588991060e-01*f[8]*z1+1.5909902576697319e+00*f[30]*z3*z4*z1+1.7787811838447134e+00*(z3*z3)*f[39]*z1+5.3033008588991060e-01*z3*f[11]+9.1855865354369182e-01*z3*f[19]*z1+9.1855865354369182e-01*z0*z1*f[16]+-5.9292706128157113e-01*f[39]*z1+-1.0269797953221864e+00*f[45]*z0*z4+1.0269797953221864e+00*(z3*z3)*f[35]+9.1855865354369182e-01*z3*z4*f[25]+1.5909902576697319e+00*f[26]*z3*z0*z1+1.7787811838447134e+00*(z3*z3)*f[42]*z4+9.1855865354369182e-01*f[21]*z0*z4+3.0618621784789724e-01*f[3]+3.0809393859665595e+00*(z3*z3)*z0*z1*f[43]+5.3033008588991060e-01*f[7]*z0+3.0809393859665595e+00*(z3*z3)*z4*z1*f[46]+-5.9292706128157113e-01*z0*f[38]+9.1855865354369182e-01*z3*f[18]*z0+-3.4232659844072882e-01*f[35]+5.3363435515341404e+00*(z3*z3)*z0*f[47]*z4*z1;

  if (dir == 3)
    return  3.5575623676894268e+00*z3*f[42]*z4*z2+2.0539595906443728e+00*z3*f[34]*z1+2.0539595906443728e+00*z3*z0*f[33]+6.1618787719331189e+00*z3*z4*z1*f[46]*z2+1.5909902576697319e+00*z0*z4*f[29]*z2+2.0539595906443728e+00*z3*z2*f[35]+5.3033008588991060e-01*z2*f[11]+3.0618621784789724e-01*f[4]+6.1618787719331189e+00*z3*z0*z1*z2*f[43]+9.1855865354369182e-01*z0*f[23]*z4+1.5909902576697319e+00*z0*z4*z1*f[28]+3.5575623676894268e+00*f[40]*z3*z0*z4+9.1855865354369182e-01*f[18]*z0*z2+1.5909902576697319e+00*f[30]*z4*z1*z2+2.0539595906443728e+00*z3*z4*f[36]+3.5575623676894268e+00*z3*z0*f[38]*z2+6.1618787719331189e+00*f[45]*z3*z0*z4*z2+3.5575623676894268e+00*z3*f[37]*z0*z1+1.0672687103068281e+01*z3*z0*f[47]*z4*z1*z2+9.1855865354369182e-01*f[19]*z1*z2+5.3033008588991060e-01*f[10]*z1+1.1858541225631423e+00*f[32]*z3+9.1855865354369182e-01*f[17]*z0*z1+3.5575623676894268e+00*z3*z4*z1*f[41]+3.5575623676894268e+00*z3*f[39]*z1*z2+5.3033008588991060e-01*z0*f[9]+9.1855865354369182e-01*z4*z2*f[25]+9.1855865354369182e-01*z4*z1*f[24]+2.7556759606310752e+00*z0*f[31]*z4*z1*z2+5.3033008588991060e-01*z4*f[15]+1.5909902576697319e+00*f[26]*z0*z1*z2+6.1618787719331189e+00*z3*z0*z4*f[44]*z1;

  if (dir == 4)
    return  3.0809393859665595e+00*(z3*z3)*z1*f[46]*z2+9.1855865354369182e-01*f[21]*z0*z2+5.3033008588991060e-01*z3*f[15]+3.0618621784789724e-01*f[5]+9.1855865354369182e-01*z3*z2*f[25]+9.1855865354369182e-01*z3*z1*f[24]+-1.0269797953221864e+00*f[45]*z0*z2+9.1855865354369182e-01*z0*z1*f[20]+1.0269797953221864e+00*(z3*z3)*f[36]+3.0809393859665595e+00*(z3*z3)*z0*f[44]*z1+-5.9292706128157113e-01*z1*f[41]+5.3363435515341404e+00*(z3*z3)*z0*f[47]*z1*z2+1.7787811838447134e+00*(z3*z3)*f[42]*z2+1.7787811838447134e+00*(z3*z3)*z1*f[41]+5.3033008588991060e-01*z0*f[12]+-5.9292706128157113e-01*f[42]*z2+2.7556759606310752e+00*z3*z0*f[31]*z1*z2+-3.4232659844072882e-01*f[36]+9.1855865354369182e-01*z3*z0*f[23]+1.5909902576697319e+00*z3*z0*f[29]*z2+9.1855865354369182e-01*f[22]*z1*z2+1.5909902576697319e+00*z0*f[27]*z1*z2+5.3033008588991060e-01*z1*f[13]+-1.0269797953221864e+00*z0*f[44]*z1+-5.9292706128157113e-01*f[40]*z0+3.0809393859665595e+00*f[45]*(z3*z3)*z0*z2+-1.7787811838447134e+00*z0*f[47]*z1*z2+5.3033008588991060e-01*f[14]*z2+1.7787811838447134e+00*f[40]*(z3*z3)*z0+1.5909902576697319e+00*z3*z0*z1*f[28]+-1.0269797953221864e+00*z1*f[46]*z2+1.5909902576697319e+00*f[30]*z3*z1*z2;

  return 0.0; // can't happen, suppresses warning

}

