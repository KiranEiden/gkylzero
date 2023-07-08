GKYL_CU_DH static inline void 
ser_3x_p2_exp_sq(const double *A, double *ASq) 
{ 
  // A:   Input DG field. 
  // ASq: Output DG field (expansion of A^2). 
 
  const double A0R2 = pow(A[0],2);
  const double A1R2 = pow(A[1],2);
  const double A2R2 = pow(A[2],2);
  const double A3R2 = pow(A[3],2);
  const double A4R2 = pow(A[4],2);
  const double A5R2 = pow(A[5],2);
  const double A6R2 = pow(A[6],2);
  const double A7R2 = pow(A[7],2);
  const double A8R2 = pow(A[8],2);
  const double A9R2 = pow(A[9],2);
  const double A10R2 = pow(A[10],2);
  const double A11R2 = pow(A[11],2);
  const double A12R2 = pow(A[12],2);
  const double A13R2 = pow(A[13],2);
  const double A14R2 = pow(A[14],2);
  const double A15R2 = pow(A[15],2);
  const double A16R2 = pow(A[16],2);
  const double A17R2 = pow(A[17],2);
  const double A18R2 = pow(A[18],2);
  const double A19R2 = pow(A[19],2);

  ASq[0] = 0.3535533905932737*A19R2+0.3535533905932737*A18R2+0.3535533905932737*A17R2+0.3535533905932737*A16R2+0.3535533905932737*A15R2+0.3535533905932737*A14R2+0.3535533905932737*A13R2+0.3535533905932737*A12R2+0.3535533905932737*A11R2+0.3535533905932737*A10R2+0.3535533905932737*A9R2+0.3535533905932737*A8R2+0.3535533905932737*A7R2+0.3535533905932737*A6R2+0.3535533905932737*A5R2+0.3535533905932737*A4R2+0.3535533905932737*A3R2+0.3535533905932737*A2R2+0.3535533905932737*A1R2+0.3535533905932737*A0R2; 
  ASq[1] = 0.7071067811865475*A[16]*A[19]+0.7071067811865475*A[14]*A[18]+0.6324555320336759*A[10]*A[17]+0.7071067811865475*A[9]*A[15]+0.632455532033676*A[5]*A[13]+0.7071067811865475*A[8]*A[12]+0.632455532033676*A[4]*A[11]+0.7071067811865475*A[6]*A[10]+0.6324555320336759*A[1]*A[7]+0.7071067811865475*A[3]*A[5]+0.7071067811865475*A[2]*A[4]+0.7071067811865475*A[0]*A[1]; 
  ASq[2] = 0.7071067811865475*A[15]*A[19]+0.6324555320336759*A[10]*A[18]+0.7071067811865475*A[13]*A[17]+0.7071067811865475*A[9]*A[16]+0.632455532033676*A[6]*A[14]+0.632455532033676*A[4]*A[12]+0.7071067811865475*A[7]*A[11]+0.7071067811865475*A[5]*A[10]+0.6324555320336759*A[2]*A[8]+0.7071067811865475*A[3]*A[6]+0.7071067811865475*A[1]*A[4]+0.7071067811865475*A[0]*A[2]; 
  ASq[3] = 0.6324555320336759*A[10]*A[19]+0.7071067811865475*A[12]*A[18]+0.7071067811865475*A[11]*A[17]+0.632455532033676*A[6]*A[16]+0.632455532033676*A[5]*A[15]+0.7071067811865475*A[8]*A[14]+0.7071067811865475*A[7]*A[13]+0.7071067811865475*A[4]*A[10]+0.6324555320336759*A[3]*A[9]+0.7071067811865475*A[2]*A[6]+0.7071067811865475*A[1]*A[5]+0.7071067811865475*A[0]*A[3]; 
  ASq[4] = 0.7071067811865475*A[9]*A[19]+0.5656854249492381*A[17]*A[18]+0.6324555320336759*A[6]*A[18]+0.6324555320336759*A[5]*A[17]+0.7071067811865475*A[15]*A[16]+0.632455532033676*A[10]*A[14]+0.632455532033676*A[10]*A[13]+0.5656854249492381*A[11]*A[12]+0.632455532033676*A[2]*A[12]+0.632455532033676*A[1]*A[11]+0.7071067811865475*A[3]*A[10]+0.6324555320336759*A[4]*A[8]+0.6324555320336759*A[4]*A[7]+0.7071067811865475*A[5]*A[6]+0.7071067811865475*A[0]*A[4]+0.7071067811865475*A[1]*A[2]; 
  ASq[5] = 0.5656854249492381*A[17]*A[19]+0.6324555320336759*A[6]*A[19]+0.7071067811865475*A[8]*A[18]+0.6324555320336759*A[4]*A[17]+0.632455532033676*A[10]*A[16]+0.5656854249492381*A[13]*A[15]+0.632455532033676*A[3]*A[15]+0.7071067811865475*A[12]*A[14]+0.632455532033676*A[1]*A[13]+0.632455532033676*A[10]*A[11]+0.7071067811865475*A[2]*A[10]+0.6324555320336759*A[5]*A[9]+0.6324555320336759*A[5]*A[7]+0.7071067811865475*A[4]*A[6]+0.7071067811865475*A[0]*A[5]+0.7071067811865475*A[1]*A[3]; 
  ASq[6] = 0.5656854249492381*A[18]*A[19]+0.6324555320336759*A[5]*A[19]+0.6324555320336759*A[4]*A[18]+0.7071067811865475*A[7]*A[17]+0.5656854249492381*A[14]*A[16]+0.632455532033676*A[3]*A[16]+0.632455532033676*A[10]*A[15]+0.632455532033676*A[2]*A[14]+0.7071067811865475*A[11]*A[13]+0.632455532033676*A[10]*A[12]+0.7071067811865475*A[1]*A[10]+0.6324555320336759*A[6]*A[9]+0.6324555320336759*A[6]*A[8]+0.7071067811865475*A[0]*A[6]+0.7071067811865475*A[4]*A[5]+0.7071067811865475*A[2]*A[3]; 
  ASq[7] = 0.3162277660168379*A19R2+0.3162277660168379*A18R2+0.2258769757263128*A17R2+0.7071067811865475*A[6]*A[17]+0.3162277660168379*A15R2+0.2258769757263128*A13R2+0.7071067811865475*A[3]*A[13]+0.3162277660168379*A12R2+0.2258769757263128*A11R2+0.7071067811865475*A[2]*A[11]+0.3162277660168379*A10R2+0.2258769757263128*A7R2+0.7071067811865475*A[0]*A[7]+0.3162277660168379*A5R2+0.3162277660168379*A4R2+0.3162277660168379*A1R2; 
  ASq[8] = 0.3162277660168379*A19R2+0.2258769757263128*A18R2+0.7071067811865475*A[5]*A[18]+0.3162277660168379*A17R2+0.3162277660168379*A16R2+0.2258769757263128*A14R2+0.7071067811865475*A[3]*A[14]+0.2258769757263128*A12R2+0.7071067811865475*A[1]*A[12]+0.3162277660168379*A11R2+0.3162277660168379*A10R2+0.2258769757263128*A8R2+0.7071067811865475*A[0]*A[8]+0.3162277660168379*A6R2+0.3162277660168379*A4R2+0.3162277660168379*A2R2; 
  ASq[9] = 0.2258769757263128*A19R2+0.7071067811865475*A[4]*A[19]+0.3162277660168379*A18R2+0.3162277660168379*A17R2+0.2258769757263128*A16R2+0.7071067811865475*A[2]*A[16]+0.2258769757263128*A15R2+0.7071067811865475*A[1]*A[15]+0.3162277660168379*A14R2+0.3162277660168379*A13R2+0.3162277660168379*A10R2+0.2258769757263128*A9R2+0.7071067811865475*A[0]*A[9]+0.3162277660168379*A6R2+0.3162277660168379*A5R2+0.3162277660168379*A3R2; 
  ASq[10] = 0.5656854249492382*A[14]*A[19]+0.5656854249492382*A[13]*A[19]+0.6324555320336759*A[3]*A[19]+0.5656854249492382*A[16]*A[18]+0.5656854249492382*A[11]*A[18]+0.6324555320336759*A[2]*A[18]+0.5656854249492382*A[15]*A[17]+0.5656854249492382*A[12]*A[17]+0.6324555320336759*A[1]*A[17]+0.632455532033676*A[5]*A[16]+0.632455532033676*A[6]*A[15]+0.632455532033676*A[4]*A[14]+0.632455532033676*A[4]*A[13]+0.632455532033676*A[6]*A[12]+0.632455532033676*A[5]*A[11]+0.6324555320336759*A[9]*A[10]+0.6324555320336759*A[8]*A[10]+0.6324555320336759*A[7]*A[10]+0.7071067811865475*A[0]*A[10]+0.7071067811865475*A[1]*A[6]+0.7071067811865475*A[2]*A[5]+0.7071067811865475*A[3]*A[4]; 
  ASq[11] = 0.6324555320336759*A[15]*A[19]+0.5656854249492382*A[10]*A[18]+0.6324555320336759*A[14]*A[17]+0.4517539514526256*A[13]*A[17]+0.7071067811865475*A[3]*A[17]+0.7071067811865475*A[6]*A[13]+0.5656854249492381*A[4]*A[12]+0.6324555320336759*A[8]*A[11]+0.4517539514526256*A[7]*A[11]+0.7071067811865475*A[0]*A[11]+0.632455532033676*A[5]*A[10]+0.7071067811865475*A[2]*A[7]+0.632455532033676*A[1]*A[4]; 
  ASq[12] = 0.6324555320336759*A[16]*A[19]+0.4517539514526256*A[14]*A[18]+0.6324555320336759*A[13]*A[18]+0.7071067811865475*A[3]*A[18]+0.5656854249492382*A[10]*A[17]+0.7071067811865475*A[5]*A[14]+0.4517539514526256*A[8]*A[12]+0.6324555320336759*A[7]*A[12]+0.7071067811865475*A[0]*A[12]+0.5656854249492381*A[4]*A[11]+0.632455532033676*A[6]*A[10]+0.7071067811865475*A[1]*A[8]+0.632455532033676*A[2]*A[4]; 
  ASq[13] = 0.5656854249492382*A[10]*A[19]+0.6324555320336759*A[12]*A[18]+0.6324555320336759*A[16]*A[17]+0.4517539514526256*A[11]*A[17]+0.7071067811865475*A[2]*A[17]+0.5656854249492381*A[5]*A[15]+0.6324555320336759*A[9]*A[13]+0.4517539514526256*A[7]*A[13]+0.7071067811865475*A[0]*A[13]+0.7071067811865475*A[6]*A[11]+0.632455532033676*A[4]*A[10]+0.7071067811865475*A[3]*A[7]+0.632455532033676*A[1]*A[5]; 
  ASq[14] = 0.5656854249492382*A[10]*A[19]+0.6324555320336759*A[15]*A[18]+0.4517539514526256*A[12]*A[18]+0.7071067811865475*A[1]*A[18]+0.6324555320336759*A[11]*A[17]+0.5656854249492381*A[6]*A[16]+0.6324555320336759*A[9]*A[14]+0.4517539514526256*A[8]*A[14]+0.7071067811865475*A[0]*A[14]+0.7071067811865475*A[5]*A[12]+0.632455532033676*A[4]*A[10]+0.7071067811865475*A[3]*A[8]+0.632455532033676*A[2]*A[6]; 
  ASq[15] = 0.4517539514526256*A[16]*A[19]+0.6324555320336759*A[11]*A[19]+0.7071067811865475*A[2]*A[19]+0.6324555320336759*A[14]*A[18]+0.5656854249492382*A[10]*A[17]+0.7071067811865475*A[4]*A[16]+0.4517539514526256*A[9]*A[15]+0.6324555320336759*A[7]*A[15]+0.7071067811865475*A[0]*A[15]+0.5656854249492381*A[5]*A[13]+0.632455532033676*A[6]*A[10]+0.7071067811865475*A[1]*A[9]+0.632455532033676*A[3]*A[5]; 
  ASq[16] = 0.4517539514526256*A[15]*A[19]+0.6324555320336759*A[12]*A[19]+0.7071067811865475*A[1]*A[19]+0.5656854249492382*A[10]*A[18]+0.6324555320336759*A[13]*A[17]+0.4517539514526256*A[9]*A[16]+0.6324555320336759*A[8]*A[16]+0.7071067811865475*A[0]*A[16]+0.7071067811865475*A[4]*A[15]+0.5656854249492381*A[6]*A[14]+0.632455532033676*A[5]*A[10]+0.7071067811865475*A[2]*A[9]+0.632455532033676*A[3]*A[6]; 
  ASq[17] = 0.5059644256269408*A[18]*A[19]+0.5656854249492381*A[5]*A[19]+0.5656854249492381*A[4]*A[18]+0.6324555320336759*A[9]*A[17]+0.6324555320336759*A[8]*A[17]+0.4517539514526256*A[7]*A[17]+0.7071067811865475*A[0]*A[17]+0.6324555320336759*A[13]*A[16]+0.5656854249492382*A[10]*A[15]+0.6324555320336759*A[11]*A[14]+0.4517539514526256*A[11]*A[13]+0.7071067811865475*A[2]*A[13]+0.5656854249492382*A[10]*A[12]+0.7071067811865475*A[3]*A[11]+0.6324555320336759*A[1]*A[10]+0.7071067811865475*A[6]*A[7]+0.6324555320336759*A[4]*A[5]; 
  ASq[18] = 0.5059644256269408*A[17]*A[19]+0.5656854249492381*A[6]*A[19]+0.6324555320336759*A[9]*A[18]+0.4517539514526256*A[8]*A[18]+0.6324555320336759*A[7]*A[18]+0.7071067811865475*A[0]*A[18]+0.5656854249492381*A[4]*A[17]+0.5656854249492382*A[10]*A[16]+0.6324555320336759*A[14]*A[15]+0.4517539514526256*A[12]*A[14]+0.7071067811865475*A[1]*A[14]+0.6324555320336759*A[12]*A[13]+0.7071067811865475*A[3]*A[12]+0.5656854249492382*A[10]*A[11]+0.6324555320336759*A[2]*A[10]+0.7071067811865475*A[5]*A[8]+0.6324555320336759*A[4]*A[6]; 
  ASq[19] = 0.4517539514526256*A[9]*A[19]+0.6324555320336759*A[8]*A[19]+0.6324555320336759*A[7]*A[19]+0.7071067811865475*A[0]*A[19]+0.5059644256269408*A[17]*A[18]+0.5656854249492381*A[6]*A[18]+0.5656854249492381*A[5]*A[17]+0.4517539514526256*A[15]*A[16]+0.6324555320336759*A[12]*A[16]+0.7071067811865475*A[1]*A[16]+0.6324555320336759*A[11]*A[15]+0.7071067811865475*A[2]*A[15]+0.5656854249492382*A[10]*A[14]+0.5656854249492382*A[10]*A[13]+0.6324555320336759*A[3]*A[10]+0.7071067811865475*A[4]*A[9]+0.6324555320336759*A[5]*A[6]; 
} 
 
