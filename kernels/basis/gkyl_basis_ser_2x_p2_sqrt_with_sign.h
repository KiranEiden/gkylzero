GKYL_CU_DH static inline void 
ser_2x_p2_sqrt_with_sign(const double *ASign, const double *A, double *ASqrt) 
{ 
  // ASign: Input DG field, used to get correct sign of Asqrt. 
  // A:     Input DG field. 
  // ASqrt: Output DG field (expansion of sqrt(A), with sign determined by Asign). 
 
  double AOrd[16] = {0.0}; 

  double temp = 0.0; 
  double temp_sign = 0.0; 
  temp = (-1.02111731798518*(A[7]+A[6]))+0.6846098004178088*(A[5]+A[4])+1.112333620718714*A[3]-0.7457659219616816*(A[2]+A[1])+0.5*A[0]; 
  temp_sign = (-1.02111731798518*(ASign[7]+ASign[6]))+0.6846098004178088*(ASign[5]+ASign[4])+1.112333620718714*ASign[3]-0.7457659219616816*(ASign[2]+ASign[1])+0.5*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[0] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[0] = -sqrt(temp); 
  } else { 
  AOrd[0] = sqrt(temp); 
  } 
  temp = 0.5446649474682879*A[7]-0.403142367494109*A[6]-0.3651715179178389*A[5]+0.6846098004178088*A[4]+0.4391550328268398*A[3]-0.2944322205496299*A[2]-0.7457659219616816*A[1]+0.5*A[0]; 
  temp_sign = 0.5446649474682879*ASign[7]-0.403142367494109*ASign[6]-0.3651715179178389*ASign[5]+0.6846098004178088*ASign[4]+0.4391550328268398*ASign[3]-0.2944322205496299*ASign[2]-0.7457659219616816*ASign[1]+0.5*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[1] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[1] = -sqrt(temp); 
  } else { 
  AOrd[1] = sqrt(temp); 
  } 
  temp = 0.5446649474682879*A[7]+0.403142367494109*A[6]-0.3651715179178389*A[5]+0.6846098004178088*A[4]-0.4391550328268398*A[3]+0.2944322205496299*A[2]-0.7457659219616816*A[1]+0.5*A[0]; 
  temp_sign = 0.5446649474682879*ASign[7]+0.403142367494109*ASign[6]-0.3651715179178389*ASign[5]+0.6846098004178088*ASign[4]-0.4391550328268398*ASign[3]+0.2944322205496299*ASign[2]-0.7457659219616816*ASign[1]+0.5*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[2] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[2] = -sqrt(temp); 
  } else { 
  AOrd[2] = sqrt(temp); 
  } 
  temp = (-1.02111731798518*A[7])+1.02111731798518*A[6]+0.6846098004178088*(A[5]+A[4])-1.112333620718714*A[3]+0.7457659219616816*A[2]-0.7457659219616816*A[1]+0.5*A[0]; 
  temp_sign = (-1.02111731798518*ASign[7])+1.02111731798518*ASign[6]+0.6846098004178088*(ASign[5]+ASign[4])-1.112333620718714*ASign[3]+0.7457659219616816*ASign[2]-0.7457659219616816*ASign[1]+0.5*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[3] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[3] = -sqrt(temp); 
  } else { 
  AOrd[3] = sqrt(temp); 
  } 
  temp = (-0.403142367494109*A[7])+0.5446649474682879*A[6]+0.6846098004178088*A[5]-0.3651715179178389*A[4]+0.4391550328268398*A[3]-0.7457659219616816*A[2]-0.2944322205496299*A[1]+0.5*A[0]; 
  temp_sign = (-0.403142367494109*ASign[7])+0.5446649474682879*ASign[6]+0.6846098004178088*ASign[5]-0.3651715179178389*ASign[4]+0.4391550328268398*ASign[3]-0.7457659219616816*ASign[2]-0.2944322205496299*ASign[1]+0.5*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[4] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[4] = -sqrt(temp); 
  } else { 
  AOrd[4] = sqrt(temp); 
  } 
  temp = 0.2150365218040566*(A[7]+A[6])-0.3651715179178389*(A[5]+A[4])+0.1733806649955719*A[3]-0.2944322205496299*(A[2]+A[1])+0.5*A[0]; 
  temp_sign = 0.2150365218040566*(ASign[7]+ASign[6])-0.3651715179178389*(ASign[5]+ASign[4])+0.1733806649955719*ASign[3]-0.2944322205496299*(ASign[2]+ASign[1])+0.5*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[5] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[5] = -sqrt(temp); 
  } else { 
  AOrd[5] = sqrt(temp); 
  } 
  temp = 0.2150365218040566*A[7]-0.2150365218040566*A[6]-0.3651715179178389*(A[5]+A[4])-0.1733806649955719*A[3]+0.2944322205496299*A[2]-0.2944322205496299*A[1]+0.5*A[0]; 
  temp_sign = 0.2150365218040566*ASign[7]-0.2150365218040566*ASign[6]-0.3651715179178389*(ASign[5]+ASign[4])-0.1733806649955719*ASign[3]+0.2944322205496299*ASign[2]-0.2944322205496299*ASign[1]+0.5*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[6] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[6] = -sqrt(temp); 
  } else { 
  AOrd[6] = sqrt(temp); 
  } 
  temp = (-0.403142367494109*A[7])-0.5446649474682879*A[6]+0.6846098004178088*A[5]-0.3651715179178389*A[4]-0.4391550328268398*A[3]+0.7457659219616816*A[2]-0.2944322205496299*A[1]+0.5*A[0]; 
  temp_sign = (-0.403142367494109*ASign[7])-0.5446649474682879*ASign[6]+0.6846098004178088*ASign[5]-0.3651715179178389*ASign[4]-0.4391550328268398*ASign[3]+0.7457659219616816*ASign[2]-0.2944322205496299*ASign[1]+0.5*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[7] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[7] = -sqrt(temp); 
  } else { 
  AOrd[7] = sqrt(temp); 
  } 
  temp = 0.403142367494109*A[7]+0.5446649474682879*A[6]+0.6846098004178088*A[5]-0.3651715179178389*A[4]-0.4391550328268398*A[3]-0.7457659219616816*A[2]+0.2944322205496299*A[1]+0.5*A[0]; 
  temp_sign = 0.403142367494109*ASign[7]+0.5446649474682879*ASign[6]+0.6846098004178088*ASign[5]-0.3651715179178389*ASign[4]-0.4391550328268398*ASign[3]-0.7457659219616816*ASign[2]+0.2944322205496299*ASign[1]+0.5*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[8] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[8] = -sqrt(temp); 
  } else { 
  AOrd[8] = sqrt(temp); 
  } 
  temp = (-0.2150365218040566*A[7])+0.2150365218040566*A[6]-0.3651715179178389*(A[5]+A[4])-0.1733806649955719*A[3]-0.2944322205496299*A[2]+0.2944322205496299*A[1]+0.5*A[0]; 
  temp_sign = (-0.2150365218040566*ASign[7])+0.2150365218040566*ASign[6]-0.3651715179178389*(ASign[5]+ASign[4])-0.1733806649955719*ASign[3]-0.2944322205496299*ASign[2]+0.2944322205496299*ASign[1]+0.5*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[9] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[9] = -sqrt(temp); 
  } else { 
  AOrd[9] = sqrt(temp); 
  } 
  temp = (-0.2150365218040566*(A[7]+A[6]))-0.3651715179178389*(A[5]+A[4])+0.1733806649955719*A[3]+0.2944322205496299*(A[2]+A[1])+0.5*A[0]; 
  temp_sign = (-0.2150365218040566*(ASign[7]+ASign[6]))-0.3651715179178389*(ASign[5]+ASign[4])+0.1733806649955719*ASign[3]+0.2944322205496299*(ASign[2]+ASign[1])+0.5*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[10] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[10] = -sqrt(temp); 
  } else { 
  AOrd[10] = sqrt(temp); 
  } 
  temp = 0.403142367494109*A[7]-0.5446649474682879*A[6]+0.6846098004178088*A[5]-0.3651715179178389*A[4]+0.4391550328268398*A[3]+0.7457659219616816*A[2]+0.2944322205496299*A[1]+0.5*A[0]; 
  temp_sign = 0.403142367494109*ASign[7]-0.5446649474682879*ASign[6]+0.6846098004178088*ASign[5]-0.3651715179178389*ASign[4]+0.4391550328268398*ASign[3]+0.7457659219616816*ASign[2]+0.2944322205496299*ASign[1]+0.5*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[11] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[11] = -sqrt(temp); 
  } else { 
  AOrd[11] = sqrt(temp); 
  } 
  temp = 1.02111731798518*A[7]-1.02111731798518*A[6]+0.6846098004178088*(A[5]+A[4])-1.112333620718714*A[3]-0.7457659219616816*A[2]+0.7457659219616816*A[1]+0.5*A[0]; 
  temp_sign = 1.02111731798518*ASign[7]-1.02111731798518*ASign[6]+0.6846098004178088*(ASign[5]+ASign[4])-1.112333620718714*ASign[3]-0.7457659219616816*ASign[2]+0.7457659219616816*ASign[1]+0.5*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[12] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[12] = -sqrt(temp); 
  } else { 
  AOrd[12] = sqrt(temp); 
  } 
  temp = (-0.5446649474682879*A[7])-0.403142367494109*A[6]-0.3651715179178389*A[5]+0.6846098004178088*A[4]-0.4391550328268398*A[3]-0.2944322205496299*A[2]+0.7457659219616816*A[1]+0.5*A[0]; 
  temp_sign = (-0.5446649474682879*ASign[7])-0.403142367494109*ASign[6]-0.3651715179178389*ASign[5]+0.6846098004178088*ASign[4]-0.4391550328268398*ASign[3]-0.2944322205496299*ASign[2]+0.7457659219616816*ASign[1]+0.5*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[13] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[13] = -sqrt(temp); 
  } else { 
  AOrd[13] = sqrt(temp); 
  } 
  temp = (-0.5446649474682879*A[7])+0.403142367494109*A[6]-0.3651715179178389*A[5]+0.6846098004178088*A[4]+0.4391550328268398*A[3]+0.2944322205496299*A[2]+0.7457659219616816*A[1]+0.5*A[0]; 
  temp_sign = (-0.5446649474682879*ASign[7])+0.403142367494109*ASign[6]-0.3651715179178389*ASign[5]+0.6846098004178088*ASign[4]+0.4391550328268398*ASign[3]+0.2944322205496299*ASign[2]+0.7457659219616816*ASign[1]+0.5*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[14] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[14] = -sqrt(temp); 
  } else { 
  AOrd[14] = sqrt(temp); 
  } 
  temp = 1.02111731798518*(A[7]+A[6])+0.6846098004178088*(A[5]+A[4])+1.112333620718714*A[3]+0.7457659219616816*(A[2]+A[1])+0.5*A[0]; 
  temp_sign = 1.02111731798518*(ASign[7]+ASign[6])+0.6846098004178088*(ASign[5]+ASign[4])+1.112333620718714*ASign[3]+0.7457659219616816*(ASign[2]+ASign[1])+0.5*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[15] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[15] = -sqrt(temp); 
  } else { 
  AOrd[15] = sqrt(temp); 
  } 
  ASqrt[0] = 0.06050149664280098*AOrd[15]+0.1134259259259259*AOrd[14]+0.1134259259259259*AOrd[13]+0.06050149664280098*AOrd[12]+0.1134259259259259*AOrd[11]+0.2126466515053472*AOrd[10]+0.2126466515053472*AOrd[9]+0.1134259259259259*AOrd[8]+0.1134259259259259*AOrd[7]+0.2126466515053472*AOrd[6]+0.2126466515053472*AOrd[5]+0.1134259259259259*AOrd[4]+0.06050149664280098*AOrd[3]+0.1134259259259259*AOrd[2]+0.1134259259259259*AOrd[1]+0.06050149664280098*AOrd[0]; 
  ASqrt[1] = 0.09023990884776016*AOrd[15]+0.1691783804450112*AOrd[14]+0.1691783804450112*AOrd[13]+0.09023990884776016*AOrd[12]+0.06679249447653647*AOrd[11]+0.1252200515903254*AOrd[10]+0.1252200515903254*AOrd[9]+0.06679249447653647*AOrd[8]-0.06679249447653647*AOrd[7]-0.1252200515903254*AOrd[6]-0.1252200515903254*AOrd[5]-0.06679249447653647*AOrd[4]-0.09023990884776016*AOrd[3]-0.1691783804450112*AOrd[2]-0.1691783804450112*AOrd[1]-0.09023990884776016*AOrd[0]; 
  ASqrt[2] = 0.09023990884776016*AOrd[15]+0.06679249447653647*AOrd[14]-0.06679249447653647*AOrd[13]-0.09023990884776016*AOrd[12]+0.1691783804450112*AOrd[11]+0.1252200515903254*AOrd[10]-0.1252200515903254*AOrd[9]-0.1691783804450112*AOrd[8]+0.1691783804450112*AOrd[7]+0.1252200515903254*AOrd[6]-0.1252200515903254*AOrd[5]-0.1691783804450112*AOrd[4]+0.09023990884776016*AOrd[3]+0.06679249447653647*AOrd[2]-0.06679249447653647*AOrd[1]-0.09023990884776016*AOrd[0]; 
  ASqrt[3] = 0.1345956976391759*AOrd[15]+0.09962313244682945*AOrd[14]-0.09962313244682945*AOrd[13]-0.1345956976391759*AOrd[12]+0.09962313244682945*AOrd[11]+0.07373763569415744*AOrd[10]-0.07373763569415744*AOrd[9]-0.09962313244682945*AOrd[8]-0.09962313244682945*AOrd[7]-0.07373763569415744*AOrd[6]+0.07373763569415744*AOrd[5]+0.09962313244682945*AOrd[4]-0.1345956976391759*AOrd[3]-0.09962313244682945*AOrd[2]+0.09962313244682945*AOrd[1]+0.1345956976391759*AOrd[0]; 
  ASqrt[4] = 0.08283983508321349*AOrd[15]+0.1553050010207067*AOrd[14]+0.1553050010207067*AOrd[13]+0.08283983508321349*AOrd[12]-0.08283983508321349*AOrd[11]-0.1553050010207067*AOrd[10]-0.1553050010207067*AOrd[9]-0.08283983508321349*AOrd[8]-0.08283983508321349*AOrd[7]-0.1553050010207067*AOrd[6]-0.1553050010207067*AOrd[5]-0.08283983508321349*AOrd[4]+0.08283983508321349*AOrd[3]+0.1553050010207067*AOrd[2]+0.1553050010207067*AOrd[1]+0.08283983508321349*AOrd[0]; 
  ASqrt[5] = 0.08283983508321349*AOrd[15]-0.08283983508321349*AOrd[14]-0.08283983508321349*AOrd[13]+0.08283983508321349*AOrd[12]+0.1553050010207067*AOrd[11]-0.1553050010207067*AOrd[10]-0.1553050010207067*AOrd[9]+0.1553050010207067*AOrd[8]+0.1553050010207067*AOrd[7]-0.1553050010207067*AOrd[6]-0.1553050010207067*AOrd[5]+0.1553050010207067*AOrd[4]+0.08283983508321349*AOrd[3]-0.08283983508321349*AOrd[2]-0.08283983508321349*AOrd[1]+0.08283983508321349*AOrd[0]; 
  ASqrt[6] = 0.1235582519719727*AOrd[15]+0.09145359262597842*AOrd[14]-0.09145359262597842*AOrd[13]-0.1235582519719727*AOrd[12]-0.1235582519719727*AOrd[11]-0.09145359262597842*AOrd[10]+0.09145359262597842*AOrd[9]+0.1235582519719727*AOrd[8]-0.1235582519719727*AOrd[7]-0.09145359262597842*AOrd[6]+0.09145359262597842*AOrd[5]+0.1235582519719727*AOrd[4]+0.1235582519719727*AOrd[3]+0.09145359262597842*AOrd[2]-0.09145359262597842*AOrd[1]-0.1235582519719727*AOrd[0]; 
  ASqrt[7] = 0.1235582519719727*AOrd[15]-0.1235582519719727*AOrd[14]-0.1235582519719727*AOrd[13]+0.1235582519719727*AOrd[12]+0.09145359262597842*AOrd[11]-0.09145359262597842*AOrd[10]-0.09145359262597842*AOrd[9]+0.09145359262597842*AOrd[8]-0.09145359262597842*AOrd[7]+0.09145359262597842*AOrd[6]+0.09145359262597842*AOrd[5]-0.09145359262597842*AOrd[4]-0.1235582519719727*AOrd[3]+0.1235582519719727*AOrd[2]+0.1235582519719727*AOrd[1]-0.1235582519719727*AOrd[0]; 

} 
 
