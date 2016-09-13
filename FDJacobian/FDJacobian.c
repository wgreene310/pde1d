/* dsm.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

typedef int integer;
typedef double real;
typedef int logical;

#include <stdlib.h>

#ifdef __cplusplus
#include <algorithm>
using std::max;
using std::min;
#else
#ifndef max
#define max(a,b) ({ a > b ? a : b; })
#endif
#ifndef min
#define min(a,b) ({ a < b ? a : b; })
#endif
#endif

#include "FDJacobian.h"


/* Table of constant values */

static integer c_n1 = -1;

int dsm_(const integer *m,  const integer *n,  const integer *npairs, 
   integer *indrow, integer *indcol, integer *ngrp, integer *maxgrp, integer *
	mingrp, integer *info, integer *ipntr, integer *jpntr, integer *iwa, 
	 const integer *liwa)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, jp, ir;
    extern /* Subroutine */ int ido_(const integer *, const integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *), seq_(const integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *), slo_(const integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);
    static integer nnz;
    extern /* Subroutine */ int degr_(const integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), setr_(const integer *, 
	    const integer *, integer *, integer *, integer *, integer *, integer *);
    static integer maxclq;
    extern /* Subroutine */ int srtdat_(const integer *, const integer *, integer *, 
	    integer *, integer *, integer *);
    static integer numgrp;
    extern /* Subroutine */ int numsrt_(const integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);

/*     ********** */

/*     SUBROUTINE DSM */

/*     THE PURPOSE OF DSM IS TO DETERMINE AN OPTIMAL OR NEAR- */
/*     OPTIMAL CONSISTENT PARTITION OF THE COLUMNS OF A SPARSE */
/*     M BY N MATRIX A. */

/*     THE SPARSITY PATTERN OF THE MATRIX A IS SPECIFIED BY */
/*     THE ARRAYS INDROW AND INDCOL. ON INPUT THE INDICES */
/*     FOR THE NON-ZERO ELEMENTS OF A ARE */

/*           INDROW(K),INDCOL(K), K = 1,2,...,NPAIRS. */

/*     THE (INDROW,INDCOL) PAIRS MAY BE SPECIFIED IN ANY ORDER. */
/*     DUPLICATE INPUT PAIRS ARE PERMITTED, BUT THE SUBROUTINE */
/*     ELIMINATES THEM. */

/*     THE SUBROUTINE PARTITIONS THE COLUMNS OF A INTO GROUPS */
/*     SUCH THAT COLUMNS IN THE SAME GROUP DO NOT HAVE A */
/*     NON-ZERO IN THE SAME ROW POSITION. A PARTITION OF THE */
/*     COLUMNS OF A WITH THIS PROPERTY IS CONSISTENT WITH THE */
/*     DIRECT DETERMINATION OF A. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE DSM(M,N,NPAIRS,INDROW,INDCOL,NGRP,MAXGRP,MINGRP, */
/*                      INFO,IPNTR,JPNTR,IWA,LIWA) */

/*     WHERE */

/*       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER */
/*         OF ROWS OF A. */

/*       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER */
/*         OF COLUMNS OF A. */

/*       NPAIRS IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE */
/*         NUMBER OF (INDROW,INDCOL) PAIRS USED TO DESCRIBE THE */
/*         SPARSITY PATTERN OF A. */

/*       INDROW IS AN INTEGER ARRAY OF LENGTH NPAIRS. ON INPUT INDROW */
/*         MUST CONTAIN THE ROW INDICES OF THE NON-ZERO ELEMENTS OF A. */
/*         ON OUTPUT INDROW IS PERMUTED SO THAT THE CORRESPONDING */
/*         COLUMN INDICES ARE IN NON-DECREASING ORDER. THE COLUMN */
/*         INDICES CAN BE RECOVERED FROM THE ARRAY JPNTR. */

/*       INDCOL IS AN INTEGER ARRAY OF LENGTH NPAIRS. ON INPUT INDCOL */
/*         MUST CONTAIN THE COLUMN INDICES OF THE NON-ZERO ELEMENTS OF */
/*         A. ON OUTPUT INDCOL IS PERMUTED SO THAT THE CORRESPONDING */
/*         ROW INDICES ARE IN NON-DECREASING ORDER. THE ROW INDICES */
/*         CAN BE RECOVERED FROM THE ARRAY IPNTR. */

/*       NGRP IS AN INTEGER OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE PARTITION OF THE COLUMNS OF A. COLUMN JCOL BELONGS */
/*         TO GROUP NGRP(JCOL). */

/*       MAXGRP IS AN INTEGER OUTPUT VARIABLE WHICH SPECIFIES THE */
/*         NUMBER OF GROUPS IN THE PARTITION OF THE COLUMNS OF A. */

/*       MINGRP IS AN INTEGER OUTPUT VARIABLE WHICH SPECIFIES A LOWER */
/*         BOUND FOR THE NUMBER OF GROUPS IN ANY CONSISTENT PARTITION */
/*         OF THE COLUMNS OF A. */

/*       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS. FOR */
/*         NORMAL TERMINATION INFO = 1. IF M, N, OR NPAIRS IS NOT */
/*         POSITIVE OR LIWA IS LESS THAN MAX(M,6*N), THEN INFO = 0. */
/*         IF THE K-TH ELEMENT OF INDROW IS NOT AN INTEGER BETWEEN */
/*         1 AND M OR THE K-TH ELEMENT OF INDCOL IS NOT AN INTEGER */
/*         BETWEEN 1 AND N, THEN INFO = -K. */

/*       IPNTR IS AN INTEGER OUTPUT ARRAY OF LENGTH M + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL. */
/*         THE COLUMN INDICES FOR ROW I ARE */

/*               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1. */

/*         NOTE THAT IPNTR(M+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       JPNTR IS AN INTEGER OUTPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW. */
/*         THE ROW INDICES FOR COLUMN J ARE */

/*               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1. */

/*         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       IWA IS AN INTEGER WORK ARRAY OF LENGTH LIWA. */

/*       LIWA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN */
/*         MAX(M,6*N). */

/*     SUBPROGRAMS CALLED */

/*       MINPACK-SUPPLIED ... DEGR,IDO,NUMSRT,SEQ,SETR,SLO,SRTDAT */

/*       FORTRAN-SUPPLIED ... MAX */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     CHECK THE INPUT DATA. */

    /* Parameter adjustments */
    --ipntr;
    --jpntr;
    --ngrp;
    --indcol;
    --indrow;
    --iwa;

    /* Function Body */
    *info = 0;
/* Computing MAX */
    i__1 = *m, i__2 = *n * 6;
    if (*m < 1 || *n < 1 || *npairs < 1 || *liwa < max(i__1,i__2)) {
	return 0;
    }
    i__1 = *npairs;
    for (k = 1; k <= i__1; ++k) {
	*info = -k;
	if (indrow[k] < 1 || indrow[k] > *m || indcol[k] < 1 || indcol[k] > *
		n) {
	    return 0;
	}
/* L10: */
    }
    *info = 1;

/*     SORT THE DATA STRUCTURE BY COLUMNS. */

    srtdat_(n, npairs, &indrow[1], &indcol[1], &jpntr[1], &iwa[1]);

/*     COMPRESS THE DATA AND DETERMINE THE NUMBER OF */
/*     NON-ZERO ELEMENTS OF A. */

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwa[i__] = 0;
/* L20: */
    }
    nnz = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	k = nnz;
	i__2 = jpntr[j + 1] - 1;
	for (jp = jpntr[j]; jp <= i__2; ++jp) {
	    ir = indrow[jp];
	    if (iwa[ir] != j) {
		indrow[nnz] = ir;
		++nnz;
		iwa[ir] = j;
	    }
/* L30: */
	}
	jpntr[j] = k;
/* L40: */
    }
    jpntr[*n + 1] = nnz;

/*     EXTEND THE DATA STRUCTURE TO ROWS. */

    setr_(m, n, &indrow[1], &jpntr[1], &indcol[1], &ipntr[1], &iwa[1]);

/*     DETERMINE A LOWER BOUND FOR THE NUMBER OF GROUPS. */

    *mingrp = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	i__2 = *mingrp, i__3 = ipntr[i__ + 1] - ipntr[i__];
	*mingrp = max(i__2,i__3);
/* L50: */
    }

/*     DETERMINE THE DEGREE SEQUENCE FOR THE INTERSECTION */
/*     GRAPH OF THE COLUMNS OF A. */

    degr_(n, &indrow[1], &jpntr[1], &indcol[1], &ipntr[1], &iwa[*n * 5 + 1], &
	    iwa[*n + 1]);

/*     COLOR THE INTERSECTION GRAPH OF THE COLUMNS OF A */
/*     WITH THE SMALLEST-LAST (SL) ORDERING. */

    slo_(n, &indrow[1], &jpntr[1], &indcol[1], &ipntr[1], &iwa[*n * 5 + 1], &
	    iwa[(*n << 2) + 1], &maxclq, &iwa[1], &iwa[*n + 1], &iwa[(*n << 1)
	     + 1], &iwa[*n * 3 + 1]);
    seq_(n, &indrow[1], &jpntr[1], &indcol[1], &ipntr[1], &iwa[(*n << 2) + 1],
	     &ngrp[1], maxgrp, &iwa[*n + 1]);
    *mingrp = max(*mingrp,maxclq);

/*     EXIT IF THE SMALLEST-LAST ORDERING IS OPTIMAL. */

    if (*maxgrp == *mingrp) {
	return 0;
    }

/*     COLOR THE INTERSECTION GRAPH OF THE COLUMNS OF A */
/*     WITH THE INCIDENCE-DEGREE (ID) ORDERING. */

    ido_(m, n, &indrow[1], &jpntr[1], &indcol[1], &ipntr[1], &iwa[*n * 5 + 1],
	     &iwa[(*n << 2) + 1], &maxclq, &iwa[1], &iwa[*n + 1], &iwa[(*n << 
	    1) + 1], &iwa[*n * 3 + 1]);
    seq_(n, &indrow[1], &jpntr[1], &indcol[1], &ipntr[1], &iwa[(*n << 2) + 1],
	     &iwa[1], &numgrp, &iwa[*n + 1]);
    *mingrp = max(*mingrp,maxclq);

/*     RETAIN THE BETTER OF THE TWO ORDERINGS SO FAR. */

    if (numgrp < *maxgrp) {
	*maxgrp = numgrp;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    ngrp[j] = iwa[j];
/* L60: */
	}

/*        EXIT IF THE INCIDENCE-DEGREE ORDERING IS OPTIMAL. */

	if (*maxgrp == *mingrp) {
	    return 0;
	}
    }

/*     COLOR THE INTERSECTION GRAPH OF THE COLUMNS OF A */
/*     WITH THE LARGEST-FIRST (LF) ORDERING. */

    i__1 = *n - 1;
    numsrt_(n, &i__1, &iwa[*n * 5 + 1], &c_n1, &iwa[(*n << 2) + 1], &iwa[(*n 
	    << 1) + 1], &iwa[*n + 1]);
    seq_(n, &indrow[1], &jpntr[1], &indcol[1], &ipntr[1], &iwa[(*n << 2) + 1],
	     &iwa[1], &numgrp, &iwa[*n + 1]);

/*     RETAIN THE BEST OF THE THREE ORDERINGS AND EXIT. */

    if (numgrp < *maxgrp) {
	*maxgrp = numgrp;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    ngrp[j] = iwa[j];
/* L70: */
	}
    }
    return 0;

/*     LAST CARD OF SUBROUTINE DSM. */

} /* dsm_ */

/* Subroutine */ int degr_(const integer *n, integer *indrow, integer *jpntr, 
	integer *indcol, integer *ipntr, integer *ndeg, integer *iwa)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer ic, ip, jp, ir, jcol;

/*     ********** */

/*     SUBROUTINE DEGR */

/*     GIVEN THE SPARSITY PATTERN OF AN M BY N MATRIX A, */
/*     THIS SUBROUTINE DETERMINES THE DEGREE SEQUENCE FOR */
/*     THE INTERSECTION GRAPH OF THE COLUMNS OF A. */

/*     IN GRAPH-THEORY TERMINOLOGY, THE INTERSECTION GRAPH OF */
/*     THE COLUMNS OF A IS THE LOOPLESS GRAPH G WITH VERTICES */
/*     A(J), J = 1,2,...,N WHERE A(J) IS THE J-TH COLUMN OF A */
/*     AND WITH EDGE (A(I),A(J)) IF AND ONLY IF COLUMNS I AND J */
/*     HAVE A NON-ZERO IN THE SAME ROW POSITION. */

/*     NOTE THAT THE VALUE OF M IS NOT NEEDED BY DEGR AND IS */
/*     THEREFORE NOT PRESENT IN THE SUBROUTINE STATEMENT. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE DEGR(N,INDROW,JPNTR,INDCOL,IPNTR,NDEG,IWA) */

/*     WHERE */

/*       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER */
/*         OF COLUMNS OF A. */

/*       INDROW IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE ROW */
/*         INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       JPNTR IS AN INTEGER INPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW. */
/*         THE ROW INDICES FOR COLUMN J ARE */

/*               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1. */

/*         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       INDCOL IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE */
/*         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       IPNTR IS AN INTEGER INPUT ARRAY OF LENGTH M + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL. */
/*         THE COLUMN INDICES FOR ROW I ARE */

/*               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1. */

/*         NOTE THAT IPNTR(M+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       NDEG IS AN INTEGER OUTPUT ARRAY OF LENGTH N WHICH */
/*         SPECIFIES THE DEGREE SEQUENCE. THE DEGREE OF THE */
/*         J-TH COLUMN OF A IS NDEG(J). */

/*       IWA IS AN INTEGER WORK ARRAY OF LENGTH N. */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     INITIALIZATION BLOCK. */

    /* Parameter adjustments */
    --iwa;
    --ndeg;
    --jpntr;
    --indrow;
    --indcol;
    --ipntr;

    /* Function Body */
    i__1 = *n;
    for (jp = 1; jp <= i__1; ++jp) {
	ndeg[jp] = 0;
	iwa[jp] = 0;
/* L10: */
    }

/*     COMPUTE THE DEGREE SEQUENCE BY DETERMINING THE CONTRIBUTIONS */
/*     TO THE DEGREES FROM THE CURRENT(JCOL) COLUMN AND FURTHER */
/*     COLUMNS WHICH HAVE NOT YET BEEN CONSIDERED. */

    i__1 = *n;
    for (jcol = 2; jcol <= i__1; ++jcol) {
	iwa[jcol] = *n;

/*        DETERMINE ALL POSITIONS (IR,JCOL) WHICH CORRESPOND */
/*        TO NON-ZEROES IN THE MATRIX. */

	i__2 = jpntr[jcol + 1] - 1;
	for (jp = jpntr[jcol]; jp <= i__2; ++jp) {
	    ir = indrow[jp];

/*           FOR EACH ROW IR, DETERMINE ALL POSITIONS (IR,IC) */
/*           WHICH CORRESPOND TO NON-ZEROES IN THE MATRIX. */

	    i__3 = ipntr[ir + 1] - 1;
	    for (ip = ipntr[ir]; ip <= i__3; ++ip) {
		ic = indcol[ip];

/*              ARRAY IWA MARKS COLUMNS WHICH HAVE CONTRIBUTED TO */
/*              THE DEGREE COUNT OF COLUMN JCOL. UPDATE THE DEGREE */
/*              COUNTS OF THESE COLUMNS AS WELL AS COLUMN JCOL. */

		if (iwa[ic] < jcol) {
		    iwa[ic] = jcol;
		    ++ndeg[ic];
		    ++ndeg[jcol];
		}
/* L20: */
	    }
/* L30: */
	}
/* L40: */
    }
    return 0;

/*     LAST CARD OF SUBROUTINE DEGR. */

} /* degr_ */

/* Subroutine */ int ido_(const integer *m, const integer *n, integer *indrow, integer *
	jpntr, integer *indcol, integer *ipntr, integer *ndeg, integer *list, 
	integer *maxclq, integer *iwa1, integer *iwa2, integer *iwa3, integer 
	*iwa4)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer ic, ip, jp, ir, jcol, ncomp, maxinc, numinc, numord, 
	    maxlst, numwgt, numlst;
    extern /* Subroutine */ int numsrt_(const integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);

/*     ********** */

/*     SUBROUTINE IDO */

/*     GIVEN THE SPARSITY PATTERN OF AN M BY N MATRIX A, THIS */
/*     SUBROUTINE DETERMINES AN INCIDENCE-DEGREE ORDERING OF THE */
/*     COLUMNS OF A. */

/*     THE INCIDENCE-DEGREE ORDERING IS DEFINED FOR THE LOOPLESS */
/*     GRAPH G WITH VERTICES A(J), J = 1,2,...,N WHERE A(J) IS THE */
/*     J-TH COLUMN OF A AND WITH EDGE (A(I),A(J)) IF AND ONLY IF */
/*     COLUMNS I AND J HAVE A NON-ZERO IN THE SAME ROW POSITION. */

/*     THE INCIDENCE-DEGREE ORDERING IS DETERMINED RECURSIVELY BY */
/*     LETTING LIST(K), K = 1,...,N BE A COLUMN WITH MAXIMAL */
/*     INCIDENCE TO THE SUBGRAPH SPANNED BY THE ORDERED COLUMNS. */
/*     AMONG ALL THE COLUMNS OF MAXIMAL INCIDENCE, IDO CHOOSES A */
/*     COLUMN OF MAXIMAL DEGREE. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE IDO(M,N,INDROW,JPNTR,INDCOL,IPNTR,NDEG,LIST, */
/*                      MAXCLQ,IWA1,IWA2,IWA3,IWA4) */

/*     WHERE */

/*       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER */
/*         OF ROWS OF A. */

/*       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER */
/*         OF COLUMNS OF A. */

/*       INDROW IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE ROW */
/*         INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       JPNTR IS AN INTEGER INPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW. */
/*         THE ROW INDICES FOR COLUMN J ARE */

/*               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1. */

/*         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       INDCOL IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE */
/*         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       IPNTR IS AN INTEGER INPUT ARRAY OF LENGTH M + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL. */
/*         THE COLUMN INDICES FOR ROW I ARE */

/*               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1. */

/*         NOTE THAT IPNTR(M+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       NDEG IS AN INTEGER INPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE DEGREE SEQUENCE. THE DEGREE OF THE J-TH COLUMN */
/*         OF A IS NDEG(J). */

/*       LIST IS AN INTEGER OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE INCIDENCE-DEGREE ORDERING OF THE COLUMNS OF A. THE J-TH */
/*         COLUMN IN THIS ORDER IS LIST(J). */

/*       MAXCLQ IS AN INTEGER OUTPUT VARIABLE SET TO THE SIZE */
/*         OF THE LARGEST CLIQUE FOUND DURING THE ORDERING. */

/*       IWA1,IWA2,IWA3, AND IWA4 ARE INTEGER WORK ARRAYS OF LENGTH N. */

/*     SUBPROGRAMS CALLED */

/*       MINPACK-SUPPLIED ... NUMSRT */

/*       FORTRAN-SUPPLIED ... MAX */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     SORT THE DEGREE SEQUENCE. */

    /* Parameter adjustments */
    --ipntr;
    --iwa4;
    --iwa3;
    --iwa2;
    --list;
    --ndeg;
    --jpntr;
    --indrow;
    --indcol;

    /* Function Body */
    i__1 = *n - 1;
    numsrt_(n, &i__1, &ndeg[1], &c_n1, &iwa4[1], &iwa2[1], &iwa3[1]);

/*     INITIALIZATION BLOCK. */

/*     CREATE A DOUBLY-LINKED LIST TO ACCESS THE INCIDENCES OF THE */
/*     COLUMNS. THE POINTERS FOR THE LINKED LIST ARE AS FOLLOWS. */

/*     EACH UN-ORDERED COLUMN IC IS IN A LIST (THE INCIDENCE LIST) */
/*     OF COLUMNS WITH THE SAME INCIDENCE. */

/*     IWA1(NUMINC) IS THE FIRST COLUMN IN THE NUMINC LIST */
/*     UNLESS IWA1(NUMINC) = 0. IN THIS CASE THERE ARE */
/*     NO COLUMNS IN THE NUMINC LIST. */

/*     IWA2(IC) IS THE COLUMN BEFORE IC IN THE INCIDENCE LIST */
/*     UNLESS IWA2(IC) = 0. IN THIS CASE IC IS THE FIRST */
/*     COLUMN IN THIS INCIDENCE LIST. */

/*     IWA3(IC) IS THE COLUMN AFTER IC IN THE INCIDENCE LIST */
/*     UNLESS IWA3(IC) = 0. IN THIS CASE IC IS THE LAST */
/*     COLUMN IN THIS INCIDENCE LIST. */

/*     IF IC IS AN UN-ORDERED COLUMN, THEN LIST(IC) IS THE */
/*     INCIDENCE OF IC TO THE GRAPH INDUCED BY THE ORDERED */
/*     COLUMNS. IF JCOL IS AN ORDERED COLUMN, THEN LIST(JCOL) */
/*     IS THE INCIDENCE-DEGREE ORDER OF COLUMN JCOL. */

    maxinc = 0;
    for (jp = *n; jp >= 1; --jp) {
	ic = iwa4[jp];
	iwa1[*n - jp] = 0;
	iwa2[ic] = 0;
	iwa3[ic] = iwa1[0];
	if (iwa1[0] > 0) {
	    iwa2[iwa1[0]] = ic;
	}
	iwa1[0] = ic;
	iwa4[jp] = 0;
	list[jp] = 0;
/* L10: */
    }

/*     DETERMINE THE MAXIMAL SEARCH LENGTH FOR THE LIST */
/*     OF COLUMNS OF MAXIMAL INCIDENCE. */

    maxlst = 0;
    i__1 = *m;
    for (ir = 1; ir <= i__1; ++ir) {
/* Computing 2nd power */
	i__2 = ipntr[ir + 1] - ipntr[ir];
	maxlst += i__2 * i__2;
/* L20: */
    }
    maxlst /= *n;
    *maxclq = 0;
    numord = 1;

/*     BEGINNING OF ITERATION LOOP. */

L30:

/*        UPDATE THE SIZE OF THE LARGEST CLIQUE */
/*        FOUND DURING THE ORDERING. */

    if (maxinc == 0) {
	ncomp = 0;
    }
    ++ncomp;
    if (maxinc + 1 == ncomp) {
	*maxclq = max(*maxclq,ncomp);
    }

/*        CHOOSE A COLUMN JCOL OF MAXIMAL DEGREE AMONG THE */
/*        COLUMNS OF MAXIMAL INCIDENCE MAXINC. */

L40:
    jp = iwa1[maxinc];
    if (jp > 0) {
	goto L50;
    }
    --maxinc;
    goto L40;
L50:
    numwgt = -1;
    i__1 = maxlst;
    for (numlst = 1; numlst <= i__1; ++numlst) {
	if (ndeg[jp] > numwgt) {
	    numwgt = ndeg[jp];
	    jcol = jp;
	}
	jp = iwa3[jp];
	if (jp <= 0) {
	    goto L70;
	}
/* L60: */
    }
L70:
    list[jcol] = numord;
    ++numord;

/*        TERMINATION TEST. */

    if (numord > *n) {
	goto L100;
    }

/*        DELETE COLUMN JCOL FROM THE MAXINC LIST. */

    if (iwa2[jcol] == 0) {
	iwa1[maxinc] = iwa3[jcol];
    } else {
	iwa3[iwa2[jcol]] = iwa3[jcol];
    }
    if (iwa3[jcol] > 0) {
	iwa2[iwa3[jcol]] = iwa2[jcol];
    }

/*        FIND ALL COLUMNS ADJACENT TO COLUMN JCOL. */

    iwa4[jcol] = *n;

/*        DETERMINE ALL POSITIONS (IR,JCOL) WHICH CORRESPOND */
/*        TO NON-ZEROES IN THE MATRIX. */

    i__1 = jpntr[jcol + 1] - 1;
    for (jp = jpntr[jcol]; jp <= i__1; ++jp) {
	ir = indrow[jp];

/*           FOR EACH ROW IR, DETERMINE ALL POSITIONS (IR,IC) */
/*           WHICH CORRESPOND TO NON-ZEROES IN THE MATRIX. */

	i__2 = ipntr[ir + 1] - 1;
	for (ip = ipntr[ir]; ip <= i__2; ++ip) {
	    ic = indcol[ip];

/*              ARRAY IWA4 MARKS COLUMNS WHICH ARE ADJACENT TO */
/*              COLUMN JCOL. */

	    if (iwa4[ic] < numord) {
		iwa4[ic] = numord;

/*                 UPDATE THE POINTERS TO THE CURRENT INCIDENCE LISTS. */

		numinc = list[ic];
		++list[ic];
/* Computing MAX */
		i__3 = maxinc, i__4 = list[ic];
		maxinc = max(i__3,i__4);

/*                 DELETE COLUMN IC FROM THE NUMINC LIST. */

		if (iwa2[ic] == 0) {
		    iwa1[numinc] = iwa3[ic];
		} else {
		    iwa3[iwa2[ic]] = iwa3[ic];
		}
		if (iwa3[ic] > 0) {
		    iwa2[iwa3[ic]] = iwa2[ic];
		}

/*                 ADD COLUMN IC TO THE NUMINC+1 LIST. */

		iwa2[ic] = 0;
		iwa3[ic] = iwa1[numinc + 1];
		if (iwa1[numinc + 1] > 0) {
		    iwa2[iwa1[numinc + 1]] = ic;
		}
		iwa1[numinc + 1] = ic;
	    }
/* L80: */
	}
/* L90: */
    }

/*        END OF ITERATION LOOP. */

    goto L30;
L100:

/*     INVERT THE ARRAY LIST. */

    i__1 = *n;
    for (jcol = 1; jcol <= i__1; ++jcol) {
	iwa2[list[jcol]] = jcol;
/* L110: */
    }
    i__1 = *n;
    for (jp = 1; jp <= i__1; ++jp) {
	list[jp] = iwa2[jp];
/* L120: */
    }
    return 0;

/*     LAST CARD OF SUBROUTINE IDO. */

} /* ido_ */

/* Subroutine */ int numsrt_(const integer *n, integer *nmax, integer *num, integer 
	*mode, integer *index, integer *last, integer *next)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l, jl, ju, jinc;

/*     **********. */

/*     SUBROUTINE NUMSRT */

/*     GIVEN A SEQUENCE OF INTEGERS, THIS SUBROUTINE GROUPS */
/*     TOGETHER THOSE INDICES WITH THE SAME SEQUENCE VALUE */
/*     AND, OPTIONALLY, SORTS THE SEQUENCE INTO EITHER */
/*     ASCENDING OR DESCENDING ORDER. */

/*     THE SEQUENCE OF INTEGERS IS DEFINED BY THE ARRAY NUM, */
/*     AND IT IS ASSUMED THAT THE INTEGERS ARE EACH FROM THE SET */
/*     0,1,...,NMAX. ON OUTPUT THE INDICES K SUCH THAT NUM(K) = L */
/*     FOR ANY L = 0,1,...,NMAX CAN BE OBTAINED FROM THE ARRAYS */
/*     LAST AND NEXT AS FOLLOWS. */

/*           K = LAST(L) */
/*           WHILE (K .NE. 0) K = NEXT(K) */

/*     OPTIONALLY, THE SUBROUTINE PRODUCES AN ARRAY INDEX SO THAT */
/*     THE SEQUENCE NUM(INDEX(I)), I = 1,2,...,N IS SORTED. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE NUMSRT(N,NMAX,NUM,MODE,INDEX,LAST,NEXT) */

/*     WHERE */

/*       N IS A POSITIVE INTEGER INPUT VARIABLE. */

/*       NMAX IS A POSITIVE INTEGER INPUT VARIABLE. */

/*       NUM IS AN INPUT ARRAY OF LENGTH N WHICH CONTAINS THE */
/*         SEQUENCE OF INTEGERS TO BE GROUPED AND SORTED. IT */
/*         IS ASSUMED THAT THE INTEGERS ARE EACH FROM THE SET */
/*         0,1,...,NMAX. */

/*       MODE IS AN INTEGER INPUT VARIABLE. THE SEQUENCE NUM IS */
/*         SORTED IN ASCENDING ORDER IF MODE IS POSITIVE AND IN */
/*         DESCENDING ORDER IF MODE IS NEGATIVE. IF MODE IS 0, */
/*         NO SORTING IS DONE. */

/*       INDEX IS AN INTEGER OUTPUT ARRAY OF LENGTH N SET SO */
/*         THAT THE SEQUENCE */

/*               NUM(INDEX(I)), I = 1,2,...,N */

/*         IS SORTED ACCORDING TO THE SETTING OF MODE. IF MODE */
/*         IS 0, INDEX IS NOT REFERENCED. */

/*       LAST IS AN INTEGER OUTPUT ARRAY OF LENGTH NMAX + 1. THE */
/*         INDEX OF NUM FOR THE LAST OCCURRENCE OF L IS LAST(L) */
/*         FOR ANY L = 0,1,...,NMAX UNLESS LAST(L) = 0. IN */
/*         THIS CASE L DOES NOT APPEAR IN NUM. */

/*       NEXT IS AN INTEGER OUTPUT ARRAY OF LENGTH N. IF */
/*         NUM(K) = L, THEN THE INDEX OF NUM FOR THE PREVIOUS */
/*         OCCURRENCE OF L IS NEXT(K) FOR ANY L = 0,1,...,NMAX */
/*         UNLESS NEXT(K) = 0. IN THIS CASE THERE IS NO PREVIOUS */
/*         OCCURRENCE OF L IN NUM. */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     DETERMINE THE ARRAYS NEXT AND LAST. */

    /* Parameter adjustments */
    --next;
    --index;
    --num;

    /* Function Body */
    i__1 = *nmax;
    for (i__ = 0; i__ <= i__1; ++i__) {
	last[i__] = 0;
/* L10: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	l = num[k];
	next[k] = last[l];
	last[l] = k;
/* L20: */
    }
    if (*mode == 0) {
	return 0;
    }

/*     STORE THE POINTERS TO THE SORTED ARRAY IN INDEX. */

    i__ = 1;
    if (*mode > 0) {
	jl = 0;
	ju = *nmax;
	jinc = 1;
    } else {
	jl = *nmax;
	ju = 0;
	jinc = -1;
    }
    i__1 = ju;
    i__2 = jinc;
    for (j = jl; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
	k = last[j];
L30:
	if (k == 0) {
	    goto L40;
	}
	index[i__] = k;
	++i__;
	k = next[k];
	goto L30;
L40:
/* L50: */
	;
    }
    return 0;

/*     LAST CARD OF SUBROUTINE NUMSRT. */

} /* numsrt_ */

/* Subroutine */ int seq_(const integer *n, integer *indrow, integer *jpntr, 
	integer *indcol, integer *ipntr, integer *list, integer *ngrp, 
	integer *maxgrp, integer *iwa)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer j, ic, ip, jp, ir, jcol;

/*     ********** */

/*     SUBROUTINE SEQ */

/*     GIVEN THE SPARSITY PATTERN OF AN M BY N MATRIX A, THIS */
/*     SUBROUTINE DETERMINES A CONSISTENT PARTITION OF THE */
/*     COLUMNS OF A BY A SEQUENTIAL ALGORITHM. */

/*     A CONSISTENT PARTITION IS DEFINED IN TERMS OF THE LOOPLESS */
/*     GRAPH G WITH VERTICES A(J), J = 1,2,...,N WHERE A(J) IS THE */
/*     J-TH COLUMN OF A AND WITH EDGE (A(I),A(J)) IF AND ONLY IF */
/*     COLUMNS I AND J HAVE A NON-ZERO IN THE SAME ROW POSITION. */

/*     A PARTITION OF THE COLUMNS OF A INTO GROUPS IS CONSISTENT */
/*     IF THE COLUMNS IN ANY GROUP ARE NOT ADJACENT IN THE GRAPH G. */
/*     IN GRAPH-THEORY TERMINOLOGY, A CONSISTENT PARTITION OF THE */
/*     COLUMNS OF A CORRESPONDS TO A COLORING OF THE GRAPH G. */

/*     THE SUBROUTINE EXAMINES THE COLUMNS IN THE ORDER SPECIFIED */
/*     BY THE ARRAY LIST, AND ASSIGNS THE CURRENT COLUMN TO THE */
/*     GROUP WITH THE SMALLEST POSSIBLE NUMBER. */

/*     NOTE THAT THE VALUE OF M IS NOT NEEDED BY SEQ AND IS */
/*     THEREFORE NOT PRESENT IN THE SUBROUTINE STATEMENT. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE SEQ(N,INDROW,JPNTR,INDCOL,IPNTR,LIST,NGRP,MAXGRP, */
/*                      IWA) */

/*     WHERE */

/*       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER */
/*         OF COLUMNS OF A. */

/*       INDROW IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE ROW */
/*         INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       JPNTR IS AN INTEGER INPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW. */
/*         THE ROW INDICES FOR COLUMN J ARE */

/*               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1. */

/*         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       INDCOL IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE */
/*         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       IPNTR IS AN INTEGER INPUT ARRAY OF LENGTH M + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL. */
/*         THE COLUMN INDICES FOR ROW I ARE */

/*               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1. */

/*         NOTE THAT IPNTR(M+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       LIST IS AN INTEGER INPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE ORDER TO BE USED BY THE SEQUENTIAL ALGORITHM. */
/*         THE J-TH COLUMN IN THIS ORDER IS LIST(J). */

/*       NGRP IS AN INTEGER OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE PARTITION OF THE COLUMNS OF A. COLUMN JCOL BELONGS */
/*         TO GROUP NGRP(JCOL). */

/*       MAXGRP IS AN INTEGER OUTPUT VARIABLE WHICH SPECIFIES THE */
/*         NUMBER OF GROUPS IN THE PARTITION OF THE COLUMNS OF A. */

/*       IWA IS AN INTEGER WORK ARRAY OF LENGTH N. */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     INITIALIZATION BLOCK. */

    /* Parameter adjustments */
    --iwa;
    --ngrp;
    --list;
    --jpntr;
    --indrow;
    --indcol;
    --ipntr;

    /* Function Body */
    *maxgrp = 0;
    i__1 = *n;
    for (jp = 1; jp <= i__1; ++jp) {
	ngrp[jp] = *n;
	iwa[jp] = 0;
/* L10: */
    }

/*     BEGINNING OF ITERATION LOOP. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jcol = list[j];

/*        FIND ALL COLUMNS ADJACENT TO COLUMN JCOL. */

/*        DETERMINE ALL POSITIONS (IR,JCOL) WHICH CORRESPOND */
/*        TO NON-ZEROES IN THE MATRIX. */

	i__2 = jpntr[jcol + 1] - 1;
	for (jp = jpntr[jcol]; jp <= i__2; ++jp) {
	    ir = indrow[jp];

/*           FOR EACH ROW IR, DETERMINE ALL POSITIONS (IR,IC) */
/*           WHICH CORRESPOND TO NON-ZEROES IN THE MATRIX. */

	    i__3 = ipntr[ir + 1] - 1;
	    for (ip = ipntr[ir]; ip <= i__3; ++ip) {
		ic = indcol[ip];

/*              ARRAY IWA MARKS THE GROUP NUMBERS OF THE */
/*              COLUMNS WHICH ARE ADJACENT TO COLUMN JCOL. */

		iwa[ngrp[ic]] = j;
/* L20: */
	    }
/* L30: */
	}

/*        ASSIGN THE SMALLEST UN-MARKED GROUP NUMBER TO JCOL. */

	i__2 = *maxgrp;
	for (jp = 1; jp <= i__2; ++jp) {
	    if (iwa[jp] != j) {
		goto L50;
	    }
/* L40: */
	}
	++(*maxgrp);
L50:
	ngrp[jcol] = jp;
/* L60: */
    }

/*        END OF ITERATION LOOP. */

    return 0;

/*     LAST CARD OF SUBROUTINE SEQ. */

} /* seq_ */

/* Subroutine */ int setr_(const integer *m, const integer *n, integer *indrow, integer *
	jpntr, integer *indcol, integer *ipntr, integer *iwa)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer jp, ir, jcol;

/*     ********** */

/*     SUBROUTINE SETR */

/*     GIVEN A COLUMN-ORIENTED DEFINITION OF THE SPARSITY PATTERN */
/*     OF AN M BY N MATRIX A, THIS SUBROUTINE DETERMINES A */
/*     ROW-ORIENTED DEFINITION OF THE SPARSITY PATTERN OF A. */

/*     ON INPUT THE COLUMN-ORIENTED DEFINITION IS SPECIFIED BY */
/*     THE ARRAYS INDROW AND JPNTR. ON OUTPUT THE ROW-ORIENTED */
/*     DEFINITION IS SPECIFIED BY THE ARRAYS INDCOL AND IPNTR. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE SETR(M,N,INDROW,JPNTR,INDCOL,IPNTR,IWA) */

/*     WHERE */

/*       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER */
/*         OF ROWS OF A. */

/*       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER */
/*         OF COLUMNS OF A. */

/*       INDROW IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE ROW */
/*         INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       JPNTR IS AN INTEGER INPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW. */
/*         THE ROW INDICES FOR COLUMN J ARE */

/*               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1. */

/*         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       INDCOL IS AN INTEGER OUTPUT ARRAY WHICH CONTAINS THE */
/*         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       IPNTR IS AN INTEGER OUTPUT ARRAY OF LENGTH M + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL. */
/*         THE COLUMN INDICES FOR ROW I ARE */

/*               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1. */

/*         NOTE THAT IPNTR(1) IS SET TO 1 AND THAT IPNTR(M+1)-1 IS */
/*         THEN THE NUMBER OF NON-ZERO ELEMENTS OF THE MATRIX A. */

/*       IWA IS AN INTEGER WORK ARRAY OF LENGTH M. */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     STORE IN ARRAY IWA THE COUNTS OF NON-ZEROES IN THE ROWS. */

    /* Parameter adjustments */
    --iwa;
    --ipntr;
    --jpntr;
    --indrow;
    --indcol;

    /* Function Body */
    i__1 = *m;
    for (ir = 1; ir <= i__1; ++ir) {
	iwa[ir] = 0;
/* L10: */
    }
    i__1 = jpntr[*n + 1] - 1;
    for (jp = 1; jp <= i__1; ++jp) {
	++iwa[indrow[jp]];
/* L20: */
    }

/*     SET POINTERS TO THE START OF THE ROWS IN INDCOL. */

    ipntr[1] = 1;
    i__1 = *m;
    for (ir = 1; ir <= i__1; ++ir) {
	ipntr[ir + 1] = ipntr[ir] + iwa[ir];
	iwa[ir] = ipntr[ir];
/* L30: */
    }

/*     FILL INDCOL. */

    i__1 = *n;
    for (jcol = 1; jcol <= i__1; ++jcol) {
	i__2 = jpntr[jcol + 1] - 1;
	for (jp = jpntr[jcol]; jp <= i__2; ++jp) {
	    ir = indrow[jp];
	    indcol[iwa[ir]] = jcol;
	    ++iwa[ir];
/* L40: */
	}
/* L50: */
    }
    return 0;

/*     LAST CARD OF SUBROUTINE SETR. */

} /* setr_ */

/* Subroutine */ int slo_(const integer *n, integer *indrow, integer *jpntr, 
	integer *indcol, integer *ipntr, integer *ndeg, integer *list, 
	integer *maxclq, integer *iwa1, integer *iwa2, integer *iwa3, integer 
	*iwa4)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer ic, ip, jp, ir, jcol, mindeg, numdeg, numord;

/*     ********** */

/*     SUBROUTINE SLO */

/*     GIVEN THE SPARSITY PATTERN OF AN M BY N MATRIX A, THIS */
/*     SUBROUTINE DETERMINES THE SMALLEST-LAST ORDERING OF THE */
/*     COLUMNS OF A. */

/*     THE SMALLEST-LAST ORDERING IS DEFINED FOR THE LOOPLESS */
/*     GRAPH G WITH VERTICES A(J), J = 1,2,...,N WHERE A(J) IS THE */
/*     J-TH COLUMN OF A AND WITH EDGE (A(I),A(J)) IF AND ONLY IF */
/*     COLUMNS I AND J HAVE A NON-ZERO IN THE SAME ROW POSITION. */

/*     THE SMALLEST-LAST ORDERING IS DETERMINED RECURSIVELY BY */
/*     LETTING LIST(K), K = N,...,1 BE A COLUMN WITH LEAST DEGREE */
/*     IN THE SUBGRAPH SPANNED BY THE UN-ORDERED COLUMNS. */

/*     NOTE THAT THE VALUE OF M IS NOT NEEDED BY SLO AND IS */
/*     THEREFORE NOT PRESENT IN THE SUBROUTINE STATEMENT. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE SLO(N,INDROW,JPNTR,INDCOL,IPNTR,NDEG,LIST, */
/*                      MAXCLQ,IWA1,IWA2,IWA3,IWA4) */

/*     WHERE */

/*       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER */
/*         OF COLUMNS OF A. */

/*       INDROW IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE ROW */
/*         INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       JPNTR IS AN INTEGER INPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW. */
/*         THE ROW INDICES FOR COLUMN J ARE */

/*               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1. */

/*         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       INDCOL IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE */
/*         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       IPNTR IS AN INTEGER INPUT ARRAY OF LENGTH M + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL. */
/*         THE COLUMN INDICES FOR ROW I ARE */

/*               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1. */

/*         NOTE THAT IPNTR(M+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       NDEG IS AN INTEGER INPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE DEGREE SEQUENCE. THE DEGREE OF THE J-TH COLUMN */
/*         OF A IS NDEG(J). */

/*       LIST IS AN INTEGER OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE SMALLEST-LAST ORDERING OF THE COLUMNS OF A. THE J-TH */
/*         COLUMN IN THIS ORDER IS LIST(J). */

/*       MAXCLQ IS AN INTEGER OUTPUT VARIABLE SET TO THE SIZE */
/*         OF THE LARGEST CLIQUE FOUND DURING THE ORDERING. */

/*       IWA1,IWA2,IWA3, AND IWA4 ARE INTEGER WORK ARRAYS OF LENGTH N. */

/*     SUBPROGRAMS CALLED */

/*       FORTRAN-SUPPLIED ... MIN */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     INITIALIZATION BLOCK. */

    /* Parameter adjustments */
    --iwa4;
    --iwa3;
    --iwa2;
    --list;
    --ndeg;
    --jpntr;
    --indrow;
    --indcol;
    --ipntr;

    /* Function Body */
    mindeg = *n;
    i__1 = *n;
    for (jp = 1; jp <= i__1; ++jp) {
	iwa1[jp - 1] = 0;
	iwa4[jp] = *n;
	list[jp] = ndeg[jp];
/* Computing MIN */
	i__2 = mindeg, i__3 = ndeg[jp];
	mindeg = min(i__2,i__3);
/* L10: */
    }

/*     CREATE A DOUBLY-LINKED LIST TO ACCESS THE DEGREES OF THE */
/*     COLUMNS. THE POINTERS FOR THE LINKED LIST ARE AS FOLLOWS. */

/*     EACH UN-ORDERED COLUMN IC IS IN A LIST (THE DEGREE LIST) */
/*     OF COLUMNS WITH THE SAME DEGREE. */

/*     IWA1(NUMDEG) IS THE FIRST COLUMN IN THE NUMDEG LIST */
/*     UNLESS IWA1(NUMDEG) = 0. IN THIS CASE THERE ARE */
/*     NO COLUMNS IN THE NUMDEG LIST. */

/*     IWA2(IC) IS THE COLUMN BEFORE IC IN THE DEGREE LIST */
/*     UNLESS IWA2(IC) = 0. IN THIS CASE IC IS THE FIRST */
/*     COLUMN IN THIS DEGREE LIST. */

/*     IWA3(IC) IS THE COLUMN AFTER IC IN THE DEGREE LIST */
/*     UNLESS IWA3(IC) = 0. IN THIS CASE IC IS THE LAST */
/*     COLUMN IN THIS DEGREE LIST. */

/*     IF IC IS AN UN-ORDERED COLUMN, THEN LIST(IC) IS THE */
/*     DEGREE OF IC IN THE GRAPH INDUCED BY THE UN-ORDERED */
/*     COLUMNS. IF JCOL IS AN ORDERED COLUMN, THEN LIST(JCOL) */
/*     IS THE SMALLEST-LAST ORDER OF COLUMN JCOL. */

    i__1 = *n;
    for (jp = 1; jp <= i__1; ++jp) {
	numdeg = ndeg[jp];
	iwa2[jp] = 0;
	iwa3[jp] = iwa1[numdeg];
	if (iwa1[numdeg] > 0) {
	    iwa2[iwa1[numdeg]] = jp;
	}
	iwa1[numdeg] = jp;
/* L20: */
    }
    *maxclq = 0;
    numord = *n;

/*     BEGINNING OF ITERATION LOOP. */

L30:

/*        MARK THE SIZE OF THE LARGEST CLIQUE */
/*        FOUND DURING THE ORDERING. */

    if (mindeg + 1 == numord && *maxclq == 0) {
	*maxclq = numord;
    }

/*        CHOOSE A COLUMN JCOL OF MINIMAL DEGREE MINDEG. */

L40:
    jcol = iwa1[mindeg];
    if (jcol > 0) {
	goto L50;
    }
    ++mindeg;
    goto L40;
L50:
    list[jcol] = numord;
    --numord;

/*        TERMINATION TEST. */

    if (numord == 0) {
	goto L80;
    }

/*        DELETE COLUMN JCOL FROM THE MINDEG LIST. */

    iwa1[mindeg] = iwa3[jcol];
    if (iwa3[jcol] > 0) {
	iwa2[iwa3[jcol]] = 0;
    }

/*        FIND ALL COLUMNS ADJACENT TO COLUMN JCOL. */

    iwa4[jcol] = 0;

/*        DETERMINE ALL POSITIONS (IR,JCOL) WHICH CORRESPOND */
/*        TO NON-ZEROES IN THE MATRIX. */

    i__1 = jpntr[jcol + 1] - 1;
    for (jp = jpntr[jcol]; jp <= i__1; ++jp) {
	ir = indrow[jp];

/*           FOR EACH ROW IR, DETERMINE ALL POSITIONS (IR,IC) */
/*           WHICH CORRESPOND TO NON-ZEROES IN THE MATRIX. */

	i__2 = ipntr[ir + 1] - 1;
	for (ip = ipntr[ir]; ip <= i__2; ++ip) {
	    ic = indcol[ip];

/*              ARRAY IWA4 MARKS COLUMNS WHICH ARE ADJACENT TO */
/*              COLUMN JCOL. */

	    if (iwa4[ic] > numord) {
		iwa4[ic] = numord;

/*                 UPDATE THE POINTERS TO THE CURRENT DEGREE LISTS. */

		numdeg = list[ic];
		--list[ic];
/* Computing MIN */
		i__3 = mindeg, i__4 = list[ic];
		mindeg = min(i__3,i__4);

/*                 DELETE COLUMN IC FROM THE NUMDEG LIST. */

		if (iwa2[ic] == 0) {
		    iwa1[numdeg] = iwa3[ic];
		} else {
		    iwa3[iwa2[ic]] = iwa3[ic];
		}
		if (iwa3[ic] > 0) {
		    iwa2[iwa3[ic]] = iwa2[ic];
		}

/*                 ADD COLUMN IC TO THE NUMDEG-1 LIST. */

		iwa2[ic] = 0;
		iwa3[ic] = iwa1[numdeg - 1];
		if (iwa1[numdeg - 1] > 0) {
		    iwa2[iwa1[numdeg - 1]] = ic;
		}
		iwa1[numdeg - 1] = ic;
	    }
/* L60: */
	}
/* L70: */
    }

/*        END OF ITERATION LOOP. */

    goto L30;
L80:

/*     INVERT THE ARRAY LIST. */

    i__1 = *n;
    for (jcol = 1; jcol <= i__1; ++jcol) {
	iwa2[list[jcol]] = jcol;
/* L90: */
    }
    i__1 = *n;
    for (jp = 1; jp <= i__1; ++jp) {
	list[jp] = iwa2[jp];
/* L100: */
    }
    return 0;

/*     LAST CARD OF SUBROUTINE SLO. */

} /* slo_ */

/* Subroutine */ int srtdat_(const integer *n, const integer *nnz, integer *indrow, 
	integer *indcol, integer *jpntr, integer *iwa)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l;

/*     ********** */

/*     SUBROUTINE SRTDAT */

/*     GIVEN THE NON-ZERO ELEMENTS OF AN M BY N MATRIX A IN */
/*     ARBITRARY ORDER AS SPECIFIED BY THEIR ROW AND COLUMN */
/*     INDICES, THIS SUBROUTINE PERMUTES THESE ELEMENTS SO */
/*     THAT THEIR COLUMN INDICES ARE IN NON-DECREASING ORDER. */

/*     ON INPUT IT IS ASSUMED THAT THE ELEMENTS ARE SPECIFIED IN */

/*           INDROW(K),INDCOL(K), K = 1,...,NNZ. */

/*     ON OUTPUT THE ELEMENTS ARE PERMUTED SO THAT INDCOL IS */
/*     IN NON-DECREASING ORDER. IN ADDITION, THE ARRAY JPNTR */
/*     IS SET SO THAT THE ROW INDICES FOR COLUMN J ARE */

/*           INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1. */

/*     NOTE THAT THE VALUE OF M IS NOT NEEDED BY SRTDAT AND IS */
/*     THEREFORE NOT PRESENT IN THE SUBROUTINE STATEMENT. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE SRTDAT(N,NNZ,INDROW,INDCOL,JPNTR,IWA) */

/*     WHERE */

/*       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER */
/*         OF COLUMNS OF A. */

/*       NNZ IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER */
/*         OF NON-ZERO ELEMENTS OF A. */

/*       INDROW IS AN INTEGER ARRAY OF LENGTH NNZ. ON INPUT INDROW */
/*         MUST CONTAIN THE ROW INDICES OF THE NON-ZERO ELEMENTS OF A. */
/*         ON OUTPUT INDROW IS PERMUTED SO THAT THE CORRESPONDING */
/*         COLUMN INDICES OF INDCOL ARE IN NON-DECREASING ORDER. */

/*       INDCOL IS AN INTEGER ARRAY OF LENGTH NNZ. ON INPUT INDCOL */
/*         MUST CONTAIN THE COLUMN INDICES OF THE NON-ZERO ELEMENTS */
/*         OF A. ON OUTPUT INDCOL IS PERMUTED SO THAT THESE INDICES */
/*         ARE IN NON-DECREASING ORDER. */

/*       JPNTR IS AN INTEGER OUTPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN THE OUTPUT */
/*         INDROW. THE ROW INDICES FOR COLUMN J ARE */

/*               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1. */

/*         NOTE THAT JPNTR(1) IS SET TO 1 AND THAT JPNTR(N+1)-1 */
/*         IS THEN NNZ. */

/*       IWA IS AN INTEGER WORK ARRAY OF LENGTH N. */

/*     SUBPROGRAMS CALLED */

/*       FORTRAN-SUPPLIED ... MAX */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     STORE IN ARRAY IWA THE COUNTS OF NON-ZEROES IN THE COLUMNS. */

    /* Parameter adjustments */
    --iwa;
    --jpntr;
    --indcol;
    --indrow;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	iwa[j] = 0;
/* L10: */
    }
    i__1 = *nnz;
    for (k = 1; k <= i__1; ++k) {
	++iwa[indcol[k]];
/* L20: */
    }

/*     SET POINTERS TO THE START OF THE COLUMNS IN INDROW. */

    jpntr[1] = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jpntr[j + 1] = jpntr[j] + iwa[j];
	iwa[j] = jpntr[j];
/* L30: */
    }
    k = 1;

/*     BEGIN IN-PLACE SORT. */

L40:
    j = indcol[k];
    if (k >= jpntr[j]) {

/*           CURRENT ELEMENT IS IN POSITION. NOW EXAMINE THE */
/*           NEXT ELEMENT OR THE FIRST UN-SORTED ELEMENT IN */
/*           THE J-TH GROUP. */

/* Computing MAX */
	i__1 = k + 1, i__2 = iwa[j];
	k = max(i__1,i__2);
    } else {

/*           CURRENT ELEMENT IS NOT IN POSITION. PLACE ELEMENT */
/*           IN POSITION AND MAKE THE DISPLACED ELEMENT THE */
/*           CURRENT ELEMENT. */

	l = iwa[j];
	++iwa[j];
	i__ = indrow[k];
	indrow[k] = indrow[l];
	indcol[k] = indcol[l];
	indrow[l] = i__;
	indcol[l] = j;
    }
    if (k <= *nnz) {
	goto L40;
    }
    return 0;

/*     LAST CARD OF SUBROUTINE SRTDAT. */

} /* srtdat_ */

int fdjs_( const integer *m,  const integer *n,  const logical *col, 
  integer *ind,
	 integer *npntr, integer *ngrp, integer *numgrp, real *d__, real *
	fjacd, real *fjac)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer ip, jp, jcol, irow;

/*     ********** */

/*     SUBROUTINE FDJS */

/*     GIVEN A CONSISTENT PARTITION OF THE COLUMNS OF AN M BY N */
/*     JACOBIAN MATRIX INTO GROUPS, THIS SUBROUTINE COMPUTES */
/*     APPROXIMATIONS TO THOSE COLUMNS IN A GIVEN GROUP.  THE */
/*     APPROXIMATIONS ARE STORED INTO EITHER A COLUMN-ORIENTED */
/*     OR A ROW-ORIENTED PATTERN. */

/*     A PARTITION IS CONSISTENT IF THE COLUMNS IN ANY GROUP */
/*     DO NOT HAVE A NON-ZERO IN THE SAME ROW POSITION. */

/*     APPROXIMATIONS TO THE COLUMNS OF THE JACOBIAN MATRIX IN A */
/*     GIVEN GROUP CAN BE OBTAINED BY SPECIFYING A DIFFERENCE */
/*     PARAMETER ARRAY D WITH D(JCOL) NON-ZERO IF AND ONLY IF */
/*     JCOL IS A COLUMN IN THE GROUP, AND AN APPROXIMATION TO */
/*     JAC*D WHERE JAC DENOTES THE JACOBIAN MATRIX OF A MAPPING F. */

/*     D CAN BE DEFINED WITH THE FOLLOWING SEGMENT OF CODE. */

/*           DO 10 JCOL = 1, N */
/*              D(JCOL) = 0.0 */
/*              IF (NGRP(JCOL) .EQ. NUMGRP) D(JCOL) = ETA(JCOL) */
/*        10    CONTINUE */

/*     IN THE ABOVE CODE NUMGRP IS THE GIVEN GROUP NUMBER, */
/*     NGRP(JCOL) IS THE GROUP NUMBER OF COLUMN JCOL, AND */
/*     ETA(JCOL) IS THE DIFFERENCE PARAMETER USED TO */
/*     APPROXIMATE COLUMN JCOL OF THE JACOBIAN MATRIX. */
/*     SUITABLE VALUES FOR THE ARRAY ETA MUST BE PROVIDED. */

/*     AS MENTIONED ABOVE, AN APPROXIMATION TO JAC*D MUST */
/*     ALSO BE PROVIDED. FOR EXAMPLE, THE APPROXIMATION */

/*           F(X+D) - F(X) */

/*     CORRESPONDS TO THE FORWARD DIFFERENCE FORMULA AT X. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE FDJS(M,N,COL,IND,NPNTR,NGRP,NUMGRP,D,FJACD,FJAC) */

/*     WHERE */

/*       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER */
/*         OF ROWS OF THE JACOBIAN MATRIX. */

/*       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER */
/*         OF COLUMNS OF THE JACOBIAN MATRIX. */

/*       COL IS A LOGICAL INPUT VARIABLE. IF COL IS SET TRUE, THEN THE */
/*         JACOBIAN APPROXIMATIONS ARE STORED INTO A COLUMN-ORIENTED */
/*         PATTERN. IF COL IS SET FALSE, THEN THE JACOBIAN */
/*         APPROXIMATIONS ARE STORED INTO A ROW-ORIENTED PATTERN. */

/*       IND IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE ROW */
/*         INDICES FOR THE NON-ZEROES IN THE JACOBIAN MATRIX */
/*         IF COL IS TRUE, AND CONTAINS THE COLUMN INDICES FOR */
/*         THE NON-ZEROES IN THE JACOBIAN MATRIX IF COL IS FALSE. */

/*       NPNTR IS AN INTEGER INPUT ARRAY WHICH SPECIFIES THE */
/*         LOCATIONS OF THE ROW INDICES IN IND IF COL IS TRUE, AND */
/*         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN IND IF */
/*         COL IS FALSE. IF COL IS TRUE, THE INDICES FOR COLUMN J ARE */

/*               IND(K), K = NPNTR(J),...,NPNTR(J+1)-1. */

/*         IF COL IS FALSE, THE INDICES FOR ROW I ARE */

/*               IND(K), K = NPNTR(I),...,NPNTR(I+1)-1. */

/*         NOTE THAT NPNTR(N+1)-1 IF COL IS TRUE, OR NPNTR(M+1)-1 */
/*         IF COL IS FALSE, IS THEN THE NUMBER OF NON-ZERO ELEMENTS */
/*         OF THE JACOBIAN MATRIX. */

/*       NGRP IS AN INTEGER INPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE PARTITION OF THE COLUMNS OF THE JACOBIAN MATRIX. */
/*         COLUMN JCOL BELONGS TO GROUP NGRP(JCOL). */

/*       NUMGRP IS A POSITIVE INTEGER INPUT VARIABLE SET TO A GROUP */
/*         NUMBER IN THE PARTITION. THE COLUMNS OF THE JACOBIAN */
/*         MATRIX IN THIS GROUP ARE TO BE ESTIMATED ON THIS CALL. */

/*       D IS AN INPUT ARRAY OF LENGTH N WHICH CONTAINS THE */
/*         DIFFERENCE PARAMETER VECTOR FOR THE ESTIMATE OF */
/*         THE JACOBIAN MATRIX COLUMNS IN GROUP NUMGRP. */

/*       FJACD IS AN INPUT ARRAY OF LENGTH M WHICH CONTAINS */
/*         AN APPROXIMATION TO THE DIFFERENCE VECTOR JAC*D, */
/*         WHERE JAC DENOTES THE JACOBIAN MATRIX. */

/*       FJAC IS AN OUTPUT ARRAY OF LENGTH NNZ, WHERE NNZ IS THE */
/*         NUMBER OF ITS NON-ZERO ELEMENTS. AT EACH CALL OF FDJS, */
/*         FJAC IS UPDATED TO INCLUDE THE NON-ZERO ELEMENTS OF THE */
/*         JACOBIAN MATRIX FOR THOSE COLUMNS IN GROUP NUMGRP. FJAC */
/*         SHOULD NOT BE ALTERED BETWEEN SUCCESSIVE CALLS TO FDJS. */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     COMPUTE ESTIMATES OF JACOBIAN MATRIX COLUMNS IN GROUP */
/*     NUMGRP. THE ARRAY FJACD MUST CONTAIN AN APPROXIMATION */
/*     TO JAC*D, WHERE JAC DENOTES THE JACOBIAN MATRIX AND D */
/*     IS A DIFFERENCE PARAMETER VECTOR WITH D(JCOL) NON-ZERO */
/*     IF AND ONLY IF JCOL IS A COLUMN IN GROUP NUMGRP. */

    /* Parameter adjustments */
    --fjacd;
    --d__;
    --ngrp;
    --ind;
    --npntr;
    --fjac;

    /* Function Body */
    if (*col) {

/*        COLUMN ORIENTATION. */

	i__1 = *n;
	for (jcol = 1; jcol <= i__1; ++jcol) {
	    if (ngrp[jcol] == *numgrp) {
		i__2 = npntr[jcol + 1] - 1;
		for (jp = npntr[jcol]; jp <= i__2; ++jp) {
		    irow = ind[jp];
		    fjac[jp] = fjacd[irow] / d__[jcol];
/* L10: */
		}
	    }
/* L20: */
	}
    } else {

/*        ROW ORIENTATION. */

	i__1 = *m;
	for (irow = 1; irow <= i__1; ++irow) {
	    i__2 = npntr[irow + 1] - 1;
	    for (ip = npntr[irow]; ip <= i__2; ++ip) {
		jcol = ind[ip];
		if (ngrp[jcol] == *numgrp) {
		    fjac[ip] = fjacd[irow] / d__[jcol];
		    goto L40;
		}
/* L30: */
	    }
L40:
/* L50: */
	    ;
	}
    }
    return 0;

/*     LAST CARD OF SUBROUTINE FDJS. */

} /* fdjs_ */

