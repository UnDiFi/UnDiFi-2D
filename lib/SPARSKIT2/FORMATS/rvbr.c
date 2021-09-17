/* f2ctmp_rvbr.f -- translated by f2c (version 19950110).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__9 = 9;
static integer c__3 = 3;
static integer c__1000 = 1000;
static integer c__6 = 6;

/* Main program */ MAIN__(void)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer maxblock;
    static doublereal a[1000], b[1000];
    static integer i, n;
    static doublereal x[100], a1[1000];
    static integer ia[101], ja[1000], ib[101], jb[1000], kb[1000], na, nc, nr,
	     nx, ny, nz, ia1[1000], ja1[1000], job;
    static doublereal ans[100];
    extern doublereal rnd_(void);
    static integer iwk[201];
    static doublereal rhs[100];
    static integer ierr;
    extern /* Subroutine */ int amux_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);
    static integer nfree, kvstc[101];
    extern /* Subroutine */ int vbrmv_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *);
    static integer kvstr[101];
    extern /* Subroutine */ int gen57bl_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, integer *, doublereal *), bsrcsr_(integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *), csrvbr_(integer *, integer *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, integer *, integer *, integer *), vbrcsr_(integer *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *);
    static doublereal stencil[700]	/* was [7][100] */;
    extern /* Subroutine */ int vbrinfo_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___19 = { 0, 6, 0, 0, 0 };
    static cilist io___30 = { 0, 6, 0, 0, 0 };
    static cilist io___31 = { 0, 6, 0, 0, 0 };
    static cilist io___32 = { 0, 6, 0, 0, 0 };
    static cilist io___34 = { 0, 6, 0, 0, 0 };


/* -----------------------------------------------------------------------
 */
/*     SPARSKIT test program for Variable Block Matrix Support */
/* -----------------------------------------------------------------------
 */
/*     This program tests all three conversion routines of csrvbr. */
/*     For each conversion to VBR, the format is converted back to CSR */
/*     with vbrcsr.  The subroutines csrkvstr, csrkvstc, and kvstmerge */
/*     become tested in the process.  The subroutines vbrinfo and vbrmv */
/*     are also tested. */
/* -----------------------------------------------------------------------
 */
/* -----dimension of grid */
    nx = 4;
    ny = 2;
    nz = 1;
    nfree = 2;
/* -----generate grid problem. */
    na = nfree * nfree;
    gen57bl_(&nx, &ny, &nz, &nfree, &na, &n, a1, ja1, ia1, iwk, stencil);
/* -----convert matrix to CSR */
    bsrcsr_(&c__1, &n, &nfree, &na, a1, ja1, ia1, a, ja, ia);
    n *= nfree;
/*     call dump(1, n, .true., a, ja, ia, 6) */
/* -----generate random x vector for testing matrix-vector product */
/* Apr. 21, 1995 */
    i__1 = n;
    for (i = 1; i <= i__1; ++i) {
	x[i - 1] = rnd_();
    }
/* -----generate correct solution for matrix-vector product */
    amux_(&n, x, ans, a, ja, ia);
    for (job = 0; job <= 2; ++job) {
	s_wsle(&io___19);
	do_lio(&c__9, &c__1, "Testing job = ", 14L);
	do_lio(&c__3, &c__1, (char *)&job, (ftnlen)sizeof(integer));
	e_wsle();
	if (job == 0) {
/* -----------maximum blocksize for random block partitioning */
	    maxblock = n / 4;
/* -----------generate random block partitioning for rows */
	    nr = 1;
	    kvstr[0] = 1;
L2000:
	    ++nr;
	    kvstr[nr - 1] = kvstr[nr - 2] + (integer) (rnd_() * maxblock) + 1;
	    if (kvstr[nr - 1] < n + 1) {
		goto L2000;
	    }
	    kvstr[nr - 1] = n + 1;
	    --nr;
/* -----------generate random block partitioning for columns */
	    nc = 1;
	    kvstc[0] = 1;
L2010:
	    ++nc;
	    kvstc[nc - 1] = kvstc[nc - 2] + (integer) (rnd_() * maxblock) + 1;
	    if (kvstc[nc - 1] < n + 1) {
		goto L2010;
	    }
	    kvstc[nc - 1] = n + 1;
	    --nc;
	}
/* --------convert to VBR format-------------------------------------
----- */
	csrvbr_(&n, ia, ja, a, &nr, &nc, kvstr, kvstc, ib, jb, kb, b, &job, 
		iwk, &c__1000, &c__1000, &ierr);
/* --------convert back to CSR format--------------------------------
----- */
	vbrcsr_(ia1, ja1, a1, &nr, kvstr, kvstc, ib, jb, kb, b, &c__1000, &
		ierr);
/* --------compare original and converted CSR structures if job not 0 
*/
	s_wsle(&io___30);
	do_lio(&c__9, &c__1, "Checking conversions....", 24L);
	e_wsle();
	if (job != 0) {
	    i__1 = n;
	    for (i = 1; i <= i__1; ++i) {
		if (ia[i - 1] != ia1[i - 1]) {
		    s_wsle(&io___31);
		    do_lio(&c__9, &c__1, "csrvbr or vbrcsr conversion mismat"
			    "ch", 36L);
		    e_wsle();
		    s_stop("", 0L);
		}
	    }
	    i__1 = ia[n] - 1;
	    for (i = 1; i <= i__1; ++i) {
		if (ja[i - 1] != ja1[i - 1] || a[i - 1] != a1[i - 1]) {
		    s_wsle(&io___32);
		    do_lio(&c__9, &c__1, "csrvbr or vbrcsr conversion mismat"
			    "ch", 36L);
		    e_wsle();
		    s_stop("", 0L);
		}
	    }
	}
/* --------test vbrinfo----------------------------------------------
----- */
	vbrinfo_(&nr, &nc, kvstr, kvstc, ib, jb, kb, iwk, &c__6);
/* --------test vbrmv------------------------------------------------
----- */
	vbrmv_(&nr, &nc, ib, jb, kb, b, kvstr, kvstc, x, rhs);
/* --------compare answer with answer computed with CSR format */
	i__1 = n;
	for (i = 1; i <= i__1; ++i) {
	    if ((d__1 = ans[i - 1] - rhs[i - 1], abs(d__1)) > (d__2 = ans[i - 
		    1] * .001, abs(d__2))) {
		s_wsle(&io___34);
		do_lio(&c__9, &c__1, "VBR matrix-vector product is erroneous "
			, 39L);
		do_lio(&c__3, &c__1, (char *)&i, (ftnlen)sizeof(integer));
		e_wsle();
		s_stop("", 0L);
	    }
	}
/* --------fill CSR structure with garbage */
	i__1 = ia1[n] - 1;
	for (i = 1; i <= i__1; ++i) {
	    ja1[i - 1] = -1;
	    a1[i - 1] = -1.;
	}
	i__1 = n + 1;
	for (i = 1; i <= i__1; ++i) {
	    ia1[i - 1] = -1;
	}
/* --------fill VBR structure with garbage */
	i__1 = kb[ib[nr] - 1] - 1;
	for (i = 1; i <= i__1; ++i) {
	    b[i - 1] = -1.;
	}
	i__1 = ib[nr];
	for (i = 1; i <= i__1; ++i) {
	    jb[i - 1] = -1;
	    kb[i - 1] = -1;
	}
	i__1 = nr + 1;
	for (i = 1; i <= i__1; ++i) {
	    ib[i - 1] = -1;
	}
/* --------fill kvstr and kvstc with garbage */
	i__1 = nr + 1;
	for (i = 1; i <= i__1; ++i) {
	    kvstr[i - 1] = -1;
	}
	i__1 = nc + 1;
	for (i = 1; i <= i__1; ++i) {
	    kvstc[i - 1] = -1;
	}
/* --------fill rhs with garbage */
	i__1 = n;
	for (i = 1; i <= i__1; ++i) {
	    rhs[i - 1] = -1.;
	}
/* -----endloop on job */
    }
    s_stop("", 0L);
    return 0;
} /* MAIN__ */

/* ----------------------------------------------------------------------- */
doublereal rnd_(void)
{
    /* Initialized data */

    static integer im = 6075;
    static integer ia = 106;
    static integer ic = 1283;
    static integer jran = 1;

    /* System generated locals */
    doublereal ret_val;

    jran = (jran * ia + ic) % im;
    ret_val = (doublereal) jran / (doublereal) im;
    return ret_val;
} /* rnd_ */

/* Main program alias */ int rvbr_ () { MAIN__ (); return 0; }
