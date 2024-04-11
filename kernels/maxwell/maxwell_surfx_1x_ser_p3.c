#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH double maxwell_surfx_1x_ser_p3(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dx1 = 2.0/dx[0]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[4]; 
  const double *ezl = &ql[8]; 
  const double *bxl = &ql[12]; 
  const double *byl = &ql[16]; 
  const double *bzl = &ql[20]; 
  const double *phl = &ql[24]; 
  const double *psl = &ql[28]; 
 
  const double *exc = &qc[0]; 
  const double *eyc = &qc[4]; 
  const double *ezc = &qc[8]; 
  const double *bxc = &qc[12]; 
  const double *byc = &qc[16]; 
  const double *bzc = &qc[20]; 
  const double *phc = &qc[24]; 
  const double *psc = &qc[28]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[4]; 
  const double *ezr = &qr[8]; 
  const double *bxr = &qr[12]; 
  const double *byr = &qr[16]; 
  const double *bzr = &qr[20]; 
  const double *phr = &qr[24]; 
  const double *psr = &qr[28]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[4]; 
  double *outEz = &out[8]; 
  double *outBx = &out[12]; 
  double *outBy = &out[16]; 
  double *outBz = &out[20]; 
  double *outPh = &out[24]; 
  double *outPs = &out[28]; 
 
  double incr_l[4]; 
 
  double incr_r[4]; 
 
  incr_l[0] = (0.6614378277661477*phl[3]-0.6614378277661477*phc[3]+0.5590169943749475*(phl[2]+phc[2])+0.4330127018922193*phl[1]-0.4330127018922193*phc[1]+0.25*(phl[0]+phc[0]))*c2chi; 
  incr_l[1] = ((-1.14564392373896*phl[3])+1.14564392373896*phc[3]-0.9682458365518543*(phl[2]+phc[2])-0.75*phl[1]+0.75*phc[1]-0.4330127018922193*(phl[0]+phc[0]))*c2chi; 
  incr_l[2] = (1.479019945774904*phl[3]-1.479019945774904*phc[3]+1.25*(phl[2]+phc[2])+0.9682458365518543*phl[1]-0.9682458365518543*phc[1]+0.5590169943749475*(phl[0]+phc[0]))*c2chi; 
  incr_l[3] = ((-1.75*phl[3])+1.75*phc[3]-1.479019945774904*(phl[2]+phc[2])-1.14564392373896*phl[1]+1.14564392373896*phc[1]-0.6614378277661477*(phl[0]+phc[0]))*c2chi; 

  incr_r[0] = (0.6614378277661477*phr[3]-0.6614378277661477*phc[3]-0.5590169943749475*(phr[2]+phc[2])+0.4330127018922193*phr[1]-0.4330127018922193*phc[1]-0.25*(phr[0]+phc[0]))*c2chi; 
  incr_r[1] = (1.14564392373896*phr[3]-1.14564392373896*phc[3]-0.9682458365518543*(phr[2]+phc[2])+0.75*phr[1]-0.75*phc[1]-0.4330127018922193*(phr[0]+phc[0]))*c2chi; 
  incr_r[2] = (1.479019945774904*phr[3]-1.479019945774904*phc[3]-1.25*(phr[2]+phc[2])+0.9682458365518543*phr[1]-0.9682458365518543*phc[1]-0.5590169943749475*(phr[0]+phc[0]))*c2chi; 
  incr_r[3] = (1.75*phr[3]-1.75*phc[3]-1.479019945774904*(phr[2]+phc[2])+1.14564392373896*phr[1]-1.14564392373896*phc[1]-0.6614378277661477*(phr[0]+phc[0]))*c2chi; 

  outEx[0] += (incr_r[0]+incr_l[0])*dx1; 
  outEx[1] += (incr_r[1]+incr_l[1])*dx1; 
  outEx[2] += (incr_r[2]+incr_l[2])*dx1; 
  outEx[3] += (incr_r[3]+incr_l[3])*dx1; 

  incr_l[0] = (0.6614378277661477*bzl[3]-0.6614378277661477*bzc[3]+0.5590169943749475*(bzl[2]+bzc[2])+0.4330127018922193*bzl[1]-0.4330127018922193*bzc[1]+0.25*(bzl[0]+bzc[0]))*c2; 
  incr_l[1] = ((-1.14564392373896*bzl[3])+1.14564392373896*bzc[3]-0.9682458365518543*(bzl[2]+bzc[2])-0.75*bzl[1]+0.75*bzc[1]-0.4330127018922193*(bzl[0]+bzc[0]))*c2; 
  incr_l[2] = (1.479019945774904*bzl[3]-1.479019945774904*bzc[3]+1.25*(bzl[2]+bzc[2])+0.9682458365518543*bzl[1]-0.9682458365518543*bzc[1]+0.5590169943749475*(bzl[0]+bzc[0]))*c2; 
  incr_l[3] = ((-1.75*bzl[3])+1.75*bzc[3]-1.479019945774904*(bzl[2]+bzc[2])-1.14564392373896*bzl[1]+1.14564392373896*bzc[1]-0.6614378277661477*(bzl[0]+bzc[0]))*c2; 

  incr_r[0] = (0.6614378277661477*bzr[3]-0.6614378277661477*bzc[3]-0.5590169943749475*(bzr[2]+bzc[2])+0.4330127018922193*bzr[1]-0.4330127018922193*bzc[1]-0.25*(bzr[0]+bzc[0]))*c2; 
  incr_r[1] = (1.14564392373896*bzr[3]-1.14564392373896*bzc[3]-0.9682458365518543*(bzr[2]+bzc[2])+0.75*bzr[1]-0.75*bzc[1]-0.4330127018922193*(bzr[0]+bzc[0]))*c2; 
  incr_r[2] = (1.479019945774904*bzr[3]-1.479019945774904*bzc[3]-1.25*(bzr[2]+bzc[2])+0.9682458365518543*bzr[1]-0.9682458365518543*bzc[1]-0.5590169943749475*(bzr[0]+bzc[0]))*c2; 
  incr_r[3] = (1.75*bzr[3]-1.75*bzc[3]-1.479019945774904*(bzr[2]+bzc[2])+1.14564392373896*bzr[1]-1.14564392373896*bzc[1]-0.6614378277661477*(bzr[0]+bzc[0]))*c2; 

  outEy[0] += (incr_r[0]+incr_l[0])*dx1; 
  outEy[1] += (incr_r[1]+incr_l[1])*dx1; 
  outEy[2] += (incr_r[2]+incr_l[2])*dx1; 
  outEy[3] += (incr_r[3]+incr_l[3])*dx1; 

  incr_l[0] = ((-0.6614378277661477*byl[3])+0.6614378277661477*byc[3]-0.5590169943749475*(byl[2]+byc[2])-0.4330127018922193*byl[1]+0.4330127018922193*byc[1]-0.25*(byl[0]+byc[0]))*c2; 
  incr_l[1] = (1.14564392373896*byl[3]-1.14564392373896*byc[3]+0.9682458365518543*(byl[2]+byc[2])+0.75*byl[1]-0.75*byc[1]+0.4330127018922193*(byl[0]+byc[0]))*c2; 
  incr_l[2] = ((-1.479019945774904*byl[3])+1.479019945774904*byc[3]-1.25*(byl[2]+byc[2])-0.9682458365518543*byl[1]+0.9682458365518543*byc[1]-0.5590169943749475*(byl[0]+byc[0]))*c2; 
  incr_l[3] = (1.75*byl[3]-1.75*byc[3]+1.479019945774904*(byl[2]+byc[2])+1.14564392373896*byl[1]-1.14564392373896*byc[1]+0.6614378277661477*(byl[0]+byc[0]))*c2; 

  incr_r[0] = ((-0.6614378277661477*byr[3])+0.6614378277661477*byc[3]+0.5590169943749475*(byr[2]+byc[2])-0.4330127018922193*byr[1]+0.4330127018922193*byc[1]+0.25*(byr[0]+byc[0]))*c2; 
  incr_r[1] = ((-1.14564392373896*byr[3])+1.14564392373896*byc[3]+0.9682458365518543*(byr[2]+byc[2])-0.75*byr[1]+0.75*byc[1]+0.4330127018922193*(byr[0]+byc[0]))*c2; 
  incr_r[2] = ((-1.479019945774904*byr[3])+1.479019945774904*byc[3]+1.25*(byr[2]+byc[2])-0.9682458365518543*byr[1]+0.9682458365518543*byc[1]+0.5590169943749475*(byr[0]+byc[0]))*c2; 
  incr_r[3] = ((-1.75*byr[3])+1.75*byc[3]+1.479019945774904*(byr[2]+byc[2])-1.14564392373896*byr[1]+1.14564392373896*byc[1]+0.6614378277661477*(byr[0]+byc[0]))*c2; 

  outEz[0] += (incr_r[0]+incr_l[0])*dx1; 
  outEz[1] += (incr_r[1]+incr_l[1])*dx1; 
  outEz[2] += (incr_r[2]+incr_l[2])*dx1; 
  outEz[3] += (incr_r[3]+incr_l[3])*dx1; 

  incr_l[0] = (0.6614378277661477*psl[3]-0.6614378277661477*psc[3]+0.5590169943749475*(psl[2]+psc[2])+0.4330127018922193*psl[1]-0.4330127018922193*psc[1]+0.25*(psl[0]+psc[0]))*gamma; 
  incr_l[1] = ((-1.14564392373896*psl[3])+1.14564392373896*psc[3]-0.9682458365518543*(psl[2]+psc[2])-0.75*psl[1]+0.75*psc[1]-0.4330127018922193*(psl[0]+psc[0]))*gamma; 
  incr_l[2] = (1.479019945774904*psl[3]-1.479019945774904*psc[3]+1.25*(psl[2]+psc[2])+0.9682458365518543*psl[1]-0.9682458365518543*psc[1]+0.5590169943749475*(psl[0]+psc[0]))*gamma; 
  incr_l[3] = ((-1.75*psl[3])+1.75*psc[3]-1.479019945774904*(psl[2]+psc[2])-1.14564392373896*psl[1]+1.14564392373896*psc[1]-0.6614378277661477*(psl[0]+psc[0]))*gamma; 

  incr_r[0] = (0.6614378277661477*psr[3]-0.6614378277661477*psc[3]-0.5590169943749475*(psr[2]+psc[2])+0.4330127018922193*psr[1]-0.4330127018922193*psc[1]-0.25*(psr[0]+psc[0]))*gamma; 
  incr_r[1] = (1.14564392373896*psr[3]-1.14564392373896*psc[3]-0.9682458365518543*(psr[2]+psc[2])+0.75*psr[1]-0.75*psc[1]-0.4330127018922193*(psr[0]+psc[0]))*gamma; 
  incr_r[2] = (1.479019945774904*psr[3]-1.479019945774904*psc[3]-1.25*(psr[2]+psc[2])+0.9682458365518543*psr[1]-0.9682458365518543*psc[1]-0.5590169943749475*(psr[0]+psc[0]))*gamma; 
  incr_r[3] = (1.75*psr[3]-1.75*psc[3]-1.479019945774904*(psr[2]+psc[2])+1.14564392373896*psr[1]-1.14564392373896*psc[1]-0.6614378277661477*(psr[0]+psc[0]))*gamma; 

  outBx[0] += (incr_r[0]+incr_l[0])*dx1; 
  outBx[1] += (incr_r[1]+incr_l[1])*dx1; 
  outBx[2] += (incr_r[2]+incr_l[2])*dx1; 
  outBx[3] += (incr_r[3]+incr_l[3])*dx1; 

  incr_l[0] = (-0.6614378277661477*ezl[3])+0.6614378277661477*ezc[3]-0.5590169943749475*(ezl[2]+ezc[2])-0.4330127018922193*ezl[1]+0.4330127018922193*ezc[1]-0.25*(ezl[0]+ezc[0]); 
  incr_l[1] = 1.14564392373896*ezl[3]-1.14564392373896*ezc[3]+0.9682458365518543*(ezl[2]+ezc[2])+0.75*ezl[1]-0.75*ezc[1]+0.4330127018922193*(ezl[0]+ezc[0]); 
  incr_l[2] = (-1.479019945774904*ezl[3])+1.479019945774904*ezc[3]-1.25*(ezl[2]+ezc[2])-0.9682458365518543*ezl[1]+0.9682458365518543*ezc[1]-0.5590169943749475*(ezl[0]+ezc[0]); 
  incr_l[3] = 1.75*ezl[3]-1.75*ezc[3]+1.479019945774904*(ezl[2]+ezc[2])+1.14564392373896*ezl[1]-1.14564392373896*ezc[1]+0.6614378277661477*(ezl[0]+ezc[0]); 

  incr_r[0] = (-0.6614378277661477*ezr[3])+0.6614378277661477*ezc[3]+0.5590169943749475*(ezr[2]+ezc[2])-0.4330127018922193*ezr[1]+0.4330127018922193*ezc[1]+0.25*(ezr[0]+ezc[0]); 
  incr_r[1] = (-1.14564392373896*ezr[3])+1.14564392373896*ezc[3]+0.9682458365518543*(ezr[2]+ezc[2])-0.75*ezr[1]+0.75*ezc[1]+0.4330127018922193*(ezr[0]+ezc[0]); 
  incr_r[2] = (-1.479019945774904*ezr[3])+1.479019945774904*ezc[3]+1.25*(ezr[2]+ezc[2])-0.9682458365518543*ezr[1]+0.9682458365518543*ezc[1]+0.5590169943749475*(ezr[0]+ezc[0]); 
  incr_r[3] = (-1.75*ezr[3])+1.75*ezc[3]+1.479019945774904*(ezr[2]+ezc[2])-1.14564392373896*ezr[1]+1.14564392373896*ezc[1]+0.6614378277661477*(ezr[0]+ezc[0]); 

  outBy[0] += (incr_r[0]+incr_l[0])*dx1; 
  outBy[1] += (incr_r[1]+incr_l[1])*dx1; 
  outBy[2] += (incr_r[2]+incr_l[2])*dx1; 
  outBy[3] += (incr_r[3]+incr_l[3])*dx1; 

  incr_l[0] = 0.6614378277661477*eyl[3]-0.6614378277661477*eyc[3]+0.5590169943749475*(eyl[2]+eyc[2])+0.4330127018922193*eyl[1]-0.4330127018922193*eyc[1]+0.25*(eyl[0]+eyc[0]); 
  incr_l[1] = (-1.14564392373896*eyl[3])+1.14564392373896*eyc[3]-0.9682458365518543*(eyl[2]+eyc[2])-0.75*eyl[1]+0.75*eyc[1]-0.4330127018922193*(eyl[0]+eyc[0]); 
  incr_l[2] = 1.479019945774904*eyl[3]-1.479019945774904*eyc[3]+1.25*(eyl[2]+eyc[2])+0.9682458365518543*eyl[1]-0.9682458365518543*eyc[1]+0.5590169943749475*(eyl[0]+eyc[0]); 
  incr_l[3] = (-1.75*eyl[3])+1.75*eyc[3]-1.479019945774904*(eyl[2]+eyc[2])-1.14564392373896*eyl[1]+1.14564392373896*eyc[1]-0.6614378277661477*(eyl[0]+eyc[0]); 

  incr_r[0] = 0.6614378277661477*eyr[3]-0.6614378277661477*eyc[3]-0.5590169943749475*(eyr[2]+eyc[2])+0.4330127018922193*eyr[1]-0.4330127018922193*eyc[1]-0.25*(eyr[0]+eyc[0]); 
  incr_r[1] = 1.14564392373896*eyr[3]-1.14564392373896*eyc[3]-0.9682458365518543*(eyr[2]+eyc[2])+0.75*eyr[1]-0.75*eyc[1]-0.4330127018922193*(eyr[0]+eyc[0]); 
  incr_r[2] = 1.479019945774904*eyr[3]-1.479019945774904*eyc[3]-1.25*(eyr[2]+eyc[2])+0.9682458365518543*eyr[1]-0.9682458365518543*eyc[1]-0.5590169943749475*(eyr[0]+eyc[0]); 
  incr_r[3] = 1.75*eyr[3]-1.75*eyc[3]-1.479019945774904*(eyr[2]+eyc[2])+1.14564392373896*eyr[1]-1.14564392373896*eyc[1]-0.6614378277661477*(eyr[0]+eyc[0]); 

  outBz[0] += (incr_r[0]+incr_l[0])*dx1; 
  outBz[1] += (incr_r[1]+incr_l[1])*dx1; 
  outBz[2] += (incr_r[2]+incr_l[2])*dx1; 
  outBz[3] += (incr_r[3]+incr_l[3])*dx1; 

  incr_l[0] = (0.6614378277661477*exl[3]-0.6614378277661477*exc[3]+0.5590169943749475*(exl[2]+exc[2])+0.4330127018922193*exl[1]-0.4330127018922193*exc[1]+0.25*(exl[0]+exc[0]))*chi; 
  incr_l[1] = ((-1.14564392373896*exl[3])+1.14564392373896*exc[3]-0.9682458365518543*(exl[2]+exc[2])-0.75*exl[1]+0.75*exc[1]-0.4330127018922193*(exl[0]+exc[0]))*chi; 
  incr_l[2] = (1.479019945774904*exl[3]-1.479019945774904*exc[3]+1.25*(exl[2]+exc[2])+0.9682458365518543*exl[1]-0.9682458365518543*exc[1]+0.5590169943749475*(exl[0]+exc[0]))*chi; 
  incr_l[3] = ((-1.75*exl[3])+1.75*exc[3]-1.479019945774904*(exl[2]+exc[2])-1.14564392373896*exl[1]+1.14564392373896*exc[1]-0.6614378277661477*(exl[0]+exc[0]))*chi; 

  incr_r[0] = (0.6614378277661477*exr[3]-0.6614378277661477*exc[3]-0.5590169943749475*(exr[2]+exc[2])+0.4330127018922193*exr[1]-0.4330127018922193*exc[1]-0.25*(exr[0]+exc[0]))*chi; 
  incr_r[1] = (1.14564392373896*exr[3]-1.14564392373896*exc[3]-0.9682458365518543*(exr[2]+exc[2])+0.75*exr[1]-0.75*exc[1]-0.4330127018922193*(exr[0]+exc[0]))*chi; 
  incr_r[2] = (1.479019945774904*exr[3]-1.479019945774904*exc[3]-1.25*(exr[2]+exc[2])+0.9682458365518543*exr[1]-0.9682458365518543*exc[1]-0.5590169943749475*(exr[0]+exc[0]))*chi; 
  incr_r[3] = (1.75*exr[3]-1.75*exc[3]-1.479019945774904*(exr[2]+exc[2])+1.14564392373896*exr[1]-1.14564392373896*exc[1]-0.6614378277661477*(exr[0]+exc[0]))*chi; 

  outPh[0] += (incr_r[0]+incr_l[0])*dx1; 
  outPh[1] += (incr_r[1]+incr_l[1])*dx1; 
  outPh[2] += (incr_r[2]+incr_l[2])*dx1; 
  outPh[3] += (incr_r[3]+incr_l[3])*dx1; 

  incr_l[0] = (0.6614378277661477*bxl[3]-0.6614378277661477*bxc[3]+0.5590169943749475*(bxl[2]+bxc[2])+0.4330127018922193*bxl[1]-0.4330127018922193*bxc[1]+0.25*(bxl[0]+bxc[0]))*c2gamma; 
  incr_l[1] = ((-1.14564392373896*bxl[3])+1.14564392373896*bxc[3]-0.9682458365518543*(bxl[2]+bxc[2])-0.75*bxl[1]+0.75*bxc[1]-0.4330127018922193*(bxl[0]+bxc[0]))*c2gamma; 
  incr_l[2] = (1.479019945774904*bxl[3]-1.479019945774904*bxc[3]+1.25*(bxl[2]+bxc[2])+0.9682458365518543*bxl[1]-0.9682458365518543*bxc[1]+0.5590169943749475*(bxl[0]+bxc[0]))*c2gamma; 
  incr_l[3] = ((-1.75*bxl[3])+1.75*bxc[3]-1.479019945774904*(bxl[2]+bxc[2])-1.14564392373896*bxl[1]+1.14564392373896*bxc[1]-0.6614378277661477*(bxl[0]+bxc[0]))*c2gamma; 

  incr_r[0] = (0.6614378277661477*bxr[3]-0.6614378277661477*bxc[3]-0.5590169943749475*(bxr[2]+bxc[2])+0.4330127018922193*bxr[1]-0.4330127018922193*bxc[1]-0.25*(bxr[0]+bxc[0]))*c2gamma; 
  incr_r[1] = (1.14564392373896*bxr[3]-1.14564392373896*bxc[3]-0.9682458365518543*(bxr[2]+bxc[2])+0.75*bxr[1]-0.75*bxc[1]-0.4330127018922193*(bxr[0]+bxc[0]))*c2gamma; 
  incr_r[2] = (1.479019945774904*bxr[3]-1.479019945774904*bxc[3]-1.25*(bxr[2]+bxc[2])+0.9682458365518543*bxr[1]-0.9682458365518543*bxc[1]-0.5590169943749475*(bxr[0]+bxc[0]))*c2gamma; 
  incr_r[3] = (1.75*bxr[3]-1.75*bxc[3]-1.479019945774904*(bxr[2]+bxc[2])+1.14564392373896*bxr[1]-1.14564392373896*bxc[1]-0.6614378277661477*(bxr[0]+bxc[0]))*c2gamma; 

  outPs[0] += (incr_r[0]+incr_l[0])*dx1; 
  outPs[1] += (incr_r[1]+incr_l[1])*dx1; 
  outPs[2] += (incr_r[2]+incr_l[2])*dx1; 
  outPs[3] += (incr_r[3]+incr_l[3])*dx1; 

  return 0.;

} 
