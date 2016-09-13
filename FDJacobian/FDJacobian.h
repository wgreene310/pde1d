// Copyright (C) 2016 William H. Greene
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses/>.

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
  int *indrow, int *indcol,
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
