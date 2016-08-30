#ifndef _FDJacobian_H_
#define _FDJacobian_H_

#if 1
#ifdef __cplusplus
extern "C" {
#endif  /* __cplusplus */
#endif

  /*
  SUBROUTINE DSM(M,N,NPAIRS,INDROW,INDCOL,NGRP,MAXGRP,MINGRP,
  *               INFO,IPNTR,JPNTR,IWA,LIWA)
  INTEGER M,N,NPAIRS,MAXGRP,MINGRP,INFO,LIWA
  INTEGER INDROW(NPAIRS),INDCOL(NPAIRS),NGRP(N),
  *        IPNTR(M+1),JPNTR(N+1),IWA(LIWA)
  */
int dsm_(const int *m, const int *n, const int *npairs, 
  const int *indrow, const int *indcol,
  int *ngrp, int *maxgrp, int *mingrp, int *info, int *ipntr,
  int *jpntr, int *iwa, const int *liwa);
/*
SUBROUTINE FDJS(M,N,COL,IND,NPNTR,NGRP,NUMGRP,D,FJACD,FJAC)
INTEGER M,N,NUMGRP
INTEGER IND(*),NPNTR(*),NGRP(N)
REAL D(N),FJACD(M),FJAC(*)
LOGICAL COL
*/
int fdjs_(const int *m, const int *n, const int *col,
  int *ind, int *npntr, int *ngrp, int *numgrp, double *d,
  double *fjacd, double *fjac);

#if 1
#ifdef __cplusplus
}

#endif  /* __cplusplus */
#endif

#endif
