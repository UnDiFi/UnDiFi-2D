!> \brief \b VKIB_SCHEME
!> \par Purpose
!>
!> This routine computes the Blended (N+LDA2) system scheme on one tetrahedron;
!> \f[
!> \Phi_i^B =  \Phi_i^N + \left( \begin{array}{llll} \alpha_1 & 0 & ... & 0 \\ 0 & \alpha_2 & ... & 0 \\ 0 & 0 & ... & \alpha_m \end{array} \right) \left( \Phi_i^{LDA_2} - \Phi_i^N \right)
!> \f]
!>
!> the distribution matrix for the LDA version 2 scheme reads:
!> \f[
!> B_i^{LDA_2} =  K_i^+ N^{-1} \quad\quad \mbox{where} \quad\quad  N = \left( \sum_{j=1}^{d+1} K_j^+ \right)
!> \f]
!>
!> so that:
!>
!> \f[
!> \Phi_i^{LDA_2}  =  -K_i^+ \left( U_{+} - U_{-} \right)
!> \f]
!>
!> whereas, for the first-order-accurate N scheme:
!>
!> \f[
!> \Phi_i^{N}  =  -K_i^+ \left( U_i - U_{-} \right)
!> \f]
!>
!> \f[
!> \Phi^{B}_i = -K_i^+ \left( U_i - U_{-} \right) + \left( \begin{array}{llll} \alpha_1 & 0 & ... & 0 \\ 0 & \alpha_2 & ... & 0 \\ 0 & 0 & ... & \alpha_m \end{array} \right) K_i^+ \left( U_i - U_{+} \right)
!> \f]
!>
!> \f[
!> \alpha = 1 - \frac{|\sum_{j=1}^{d+1} \phi_j^N|}{\sum_{j=1}^{d+1} |\phi_j^N|}
!> \f]
!>
!>
!>
!> @param[in] MATRIXSPLITTER is the subroutine use to split the Jacobian matrix
!> @param[in] W is the set of NORDER dependent variables in the NOFVERT vertices of the current cell
!> @param[out] NODRES are the signals sent to each of the NOFVERT vertices
!> @param[in,out] TSTEP is the elemental contribution to the (inverse of the) time step
!> @param[out] BETA the NOFVERT distribution matrices
!> @param[out] STIFC the elemental contribution to the implicit matrix (computed only if MATRIX_ASSEMBLY is true)
!> @param[in] NORDER is the order of the hyperbolic system being solved, it could be < than NDOF
!> @param[in] NDOF is the number of dofs in the gridpoints, also the leading dimension of W, STIFC, NODRES, TSTEP
!> @param[out] NOFVERT is number of vertices of the simplicial elements (= NDIM+1)
!> @param[in] VCN the NDIM cartesian component of the inward face normal to all NOFVERT vertices, scaled by its measure
!> @param[in] NDIM is the dimension of the space
!> @param[in] DFGHDU the NDIM cartesian components of the Jacobian matrix
!> @param[in] LDJ is the leading dimension of DFGHDU
!> @param[out] CELRES the elemental residual \f$ = -\Phi = \sum_{k=j}^{d+1} K_j U_j \f$
!> @param[in] SOURCE the elemental source term \f$ = -\int S \mathrm{d}V \f$
!> @param[in] IELEM is the current simplicial element
!> @param[in] MATRIX_ASSEMBLY when set == .TRUE. the STIFC matrix will be assembled
!> \author $Author: abonfi $
!> \version $Revision: 1.12 $
!> \date $Date: 2020/03/28 09:49:28 $

      SUBROUTINE VKIB_SCHEME(MATRIXSPLITTER,W,NODRES,TSTEP,BETA,STIFC,
     +                       NORDER,NDOF,NOFVERT,VCN,NDIM,DFGHDU,LDJ,
     +                       CELRES,SOURCE,IELEM,MATRIX_ASSEMBLY)
C
      IMPLICIT NONE
C
C     $Id: VKIB_scheme.F,v 1.12 2020/03/28 09:49:28 abonfi Exp $
C
C
C
C
C This routine computes the system VKI-B scheme on one tetrahedron
C
      include 'paramt.h'
      include 'constants.h'
      include 'flags.com'
      include 'three.com'
C
C When DEBUG is defined, you should also use the runtime option -check
C
CCC#define DEBUG
C
C When DUMP_WEIGHT is defined, the weights are written in fort.12 to be later examined using celplot
C When using DUMP_WEIGHT, you should first converge the calculation and then restart from the solution and run only a pseudo-time-step
C
CCCC#define DUMP_WEIGHT
C
C These two defines are used to check that the B scheme reverts to either LDA or N, when ALFA = 1, resp. 0
C
CC#define N_SCHEME
CC#define LDA_SCHEME
C
C     .. Parameters ..
C     ..
C     .. Scalar Arguments ..
      INTEGER IELEM,LDJ,NDIM,NOFVERT,NDOF,NORDER
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DFGHDU(LDJ,*),NODRES(NDOF,NOFVERT),
     +                 TSTEP(NDOF,NOFVERT),SOURCE(*),
     +                 CELRES(*),BETA(*),
     4                 STIFC(NDOF,NDOF,NOFVERT,NOFVERT),
     +                 VCN(NDIM,NOFVERT),W(NDOF,NOFVERT)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL MATRIXSPLITTER
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION HELP,DUM,TEMP,TMP
      INTEGER I,INFO,IROW,IVAR,IVERT,J,JCOL
      LOGICAL LFLAG
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION KMAT(MAXNOFVAR,MAXNOFVAR),
     +KNEG(MAXNOFVAR,MAXNOFVAR,VMAX),KPOS(MAXNOFVAR,MAXNOFVAR,VMAX),
     +SUM_K_NEG(MAXNOFVAR,MAXNOFVAR),UNEG(MAXNOFVAR),UPOS(MAXNOFVAR),
     +VLEFT(MAXNOFVAR,MAXNOFVAR),VRIGHT(MAXNOFVAR,MAXNOFVAR),
     +WKSP1(MAXNOFVAR),WKSP2(MAXNOFVAR),WKSP3(MAXNOFVAR,MAXNOFVAR),
     +WNEG(MAXNOFVAR),WPOS(MAXNOFVAR),WR(MAXNOFVAR),
     +PHIN(MAXNOFVAR,VMAX),DPHI(MAXNOFVAR,VMAX),PHI(MAXNOFVAR),
     +ALFA(MAXNOFVAR),SUM_K_POS(MAXNOFVAR,MAXNOFVAR),WKSP4(MAXNOFVAR)
#ifdef DEBUG
      DOUBLE PRECISION DW(MAXNOFVAR*VMAX)
#endif
      INTEGER IPIV(MAXNOFVAR),JPIV(MAXNOFVAR)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DNRM2
      LOGICAL NULLMAT,UNITMAT
      EXTERNAL DNRM2,NULLMAT,UNITMAT
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DGEMM,DGEMV,DGETRF,DGETRS,DINIT,MATSUB,
     +         MATSUM,R8Mat_Print
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS
C     ..
C     .. Statement Functions ..
C     ..
      DO 1 J = 1,MAXNOFVAR
          UNEG(J) = ZERO
          UPOS(J) = ZERO
          DO 1 I = 1,MAXNOFVAR
              SUM_K_NEG(I,J) = ZERO
              SUM_K_POS(I,J) = ZERO
    1 CONTINUE
      DO 3 I = 1,2*NDOF
          CELRES(I) = ZERO
    3 CONTINUE
C
C --------------- Debugging code starts here ---------------
C
      IF (ICHECK.NE.0) THEN
          CALL DINIT(MAXNOFVAR,ZERO,WKSP1,1)
          CALL DINIT(MAXNOFVAR,ZERO,WKSP2,1)
          CALL DINIT(MAXNOFVAR**2,ZERO,WKSP3,1)
      ENDIF
C
C --------------- Debugging code ends here ---------------
C
      DO 10 IVERT = 1,NOFVERT
C
C
C    The matrix is split into its positive and negative parts
C
C    Numerically or Analitically ..
C
          CALL MATRIXSPLITTER(IELEM,NDIM,NORDER,VCN(1,IVERT),DFGHDU,LDJ,
     +         KMAT,KPOS(1,1,IVERT),KNEG(1,1,IVERT),
     +         VLEFT,VRIGHT,MAXNOFVAR,WR,WPOS,WNEG,.TRUE.)
C
C    The negative jacobians are accumulated in SUM_K_NEG ..
C
          DO 15 JCOL = 1,NORDER
              DO 15 IROW = 1,NORDER
                  SUM_K_NEG(IROW,JCOL) = SUM_K_NEG(IROW,JCOL) +
     +                                   KNEG(IROW,JCOL,IVERT)
                  SUM_K_POS(IROW,JCOL) = SUM_K_POS(IROW,JCOL) +
     +                                   KPOS(IROW,JCOL,IVERT)
   15     CONTINUE
C
C
C    Computes KNEG(IVERT)*U(IVERT) and adds it to UNEG ..
C
cblas    CALL DGEMV('N',NORDER,NORDER,ONE,KNEG(1,1,IVERT),MAXNOFVAR,W
cblas+   (1,IVERT),1,ONE,UNEG,1)
C
          DO 11 I = 1,NORDER
              HELP = ZERO
              DUM = ZERO
              DO 13 J = 1,NORDER
                  HELP = HELP + KNEG(I,J,IVERT)*W(J,IVERT)
                  DUM = DUM + KPOS(I,J,IVERT)*W(J,IVERT)
   13         CONTINUE
              UNEG(I) = UNEG(I) + HELP
C    we use a (-) here since we shall later multiply by SUM_K_NEG^{-1}
              UPOS(I) = UPOS(I) - DUM
caldo
   11     CONTINUE
C
C
C --------------- Debugging code starts here ---------------
C
          IF (ICHECK.NE.0) THEN
C
C    the residual is computes as Sum_j K(j) * U(j) ..
C
              CALL DGEMV('N',NORDER,NORDER,ONE,KMAT,MAXNOFVAR,
     +                   W(1,IVERT),1,ONE,CELRES,1)
          ENDIF
C
C --------------- Debugging code ends here ---------------
C
C    Timestep ..
C
          DO 34 I = 1,NORDER
              TSTEP(I,IVERT) = TSTEP(I,IVERT) + WPOS(I)
   34     CONTINUE
C
   10 CONTINUE ! End loop over vertices
C
C    Finds the generalized inflow points
C       solving [\sum_j K_j^-] U_{\pm} = \sum_j K_j^- U_j
C
C    LU factorization ..
C
cnag  CALL F07ADF(NORDER,NORDER,SUM_K_NEG,MAXNOFVAR,IPIV,INFO)
      CALL DGETRF(NORDER,NORDER,SUM_K_NEG,MAXNOFVAR,IPIV,INFO)
caldo
caldo CALL DGETRF(NORDER,NORDER,SUM_K_POS,MAXNOFVAR,JPIV,INFO)
caldo
C
      IF (INFO.GT.0) THEN
          WRITE (6,FMT=99999) IELEM

99999     FORMAT (5X,'Matrix SUM_K_NEG is singular in IELEM = ',I6)
*
          DO 9 IVERT = 1,NOFVERT
              WRITE (6,FMT=*) '    Vertice # ',IVERT
              CALL R8Mat_Print('G',' ',NORDER,NORDER,KPOS(1,1,IVERT),
     +                    MAXNOFVAR,'K(+) ',INFO)
              CALL R8Mat_Print('G',' ',NORDER,NORDER,KNEG(1,1,IVERT),
     +                    MAXNOFVAR,'K(-) ',INFO)
    9     CONTINUE
          CALL R8Mat_Print('G',' ',NORDER,NORDER,SUM_K_NEG,MAXNOFVAR,
     +                'SUM_j K(-) ',INFO)
          STOP

      ENDIF
C
C    solution ..
C
cnag  CALL F07AEF('N',NORDER,1,SUM_K_NEG,MAXNOFVAR,IPIV,UNEG,MAXNOFVAR,INFO)
      CALL DGETRS('N',NORDER,1,SUM_K_NEG,MAXNOFVAR,IPIV,UNEG,
     +MAXNOFVAR,INFO)
      CALL DGETRS('N',NORDER,1,SUM_K_NEG,MAXNOFVAR,IPIV,UPOS,
     +MAXNOFVAR,INFO)
caldo
#ifdef DEBUG
caldo
caldo check - sum_k_pos (u+ - u-) = \phi
caldo
caldo
caldo
      call dcopy(norder,upos,1,DW(norder+1),1)
      call daxpy(norder,-one,uneg,1,DW(norder+1),1) ! (u+ - u-)
      call dcopy(norder,CELRES,1,DW(1),1)
      call dgemv('N',norder,norder,mone,SUM_K_POS,MAXNOFVAR,
     &DW(norder+1),1,ONE,DW(1),1) ! compute the difference btw. \Phi and sum_k_pos (u+ - u-)
      write(6,*)'e = ',ielem
      write(6,*)' \Phi - sum_k_pos (u+ - u-) = ',(DW(i),i=1,norder)
      write(6,*)' \Phi = ',(celres(i),i=1,norder)
#endif
C
      DO IVAR = 1,NORDER
         PHI(IVAR) = ZERO
         ALFA(IVAR) = ZERO
      ENDDO
C
C Loops over nodes
C
      DO 30 IVERT = 1,NOFVERT
!         I = (IVERT-1)*NORDER
C
C    NODRES[INODE] = - KPOS(IVERT)*[U(IVERT)-UNEG]
C
cblas    CALL DGEMV('N',NORDER,NORDER,ONE,KPOS(1,1,IVERT),MAXNOFVAR, W
cblas+   (1,IVERT),1,ZERO,WKSP,1)
C
          DO 23 IVAR = 1,NORDER
              HELP = ZERO
              TMP = ZERO
              DO 21 J = 1,NORDER
                  TEMP = KPOS(IVAR,J,IVERT)
                  HELP = HELP + TEMP * (W(J,IVERT)-UNEG(J)) ! N-scheme contribution
                  TMP = TMP - TEMP * (W(J,IVERT)-UPOS(J)) ! anti-diffusive correction
   21         CONTINUE
              PHIN(IVAR,IVERT) = -HELP
              DPHI(IVAR,IVERT) = -TMP
               PHI(IVAR) =  PHI(IVAR) + PHIN(IVAR,IVERT)
              ALFA(IVAR) = ALFA(IVAR) + ABS(PHIN(IVAR,IVERT))
   23     CONTINUE

C
C --------------- Debugging code starts here ---------------
C
          IF (ICHECK.NE.0) THEN
C
C    The residual is computed as Sum_j KPOS(j) * [U(j)-Uin]
C       and stored in WKSP1
C

              CALL DAXPY(NORDER,-ONE,PHIN(1,IVERT),1,WKSP1,1)
          ENDIF
C
C --------------- Debugging code ends here ---------------
C
   30 CONTINUE ! End loop over nodes
C
C
      DO 32 IVAR = 1,NORDER
          IF((ALFA(IVAR)).EQ.ZERO)THEN
                  ALFA(IVAR) = ONE
          ELSE
              ALFA(IVAR) = ONE - ABS(PHI(IVAR))/ALFA(IVAR)
          ENDIF
          if( (alfa(ivar).GT.ONE) .OR. (alfa(ivar).LT.ZERO) )then
                 write(6,*)'Smthg wrong with alfa ',ivar,ALFA(IVAR)
          endif
   32 CONTINUE ! end loop over variables
#ifdef N_SCHEME
      DO IVAR = 1, NORDER
         ALFA(IVAR) = ZERO
      ENDDO
#endif
#ifdef LDA_SCHEME
      DO IVAR = 1, NORDER
         ALFA(IVAR) = ONE
      ENDDO
#endif
#ifdef DUMP_WEIGHT
      write(12,*)ielem,(ALFA(IVAR),IVAR = 1,NORDER)
#endif
C
C     compute the blended update
C
      IF(ICHECK.NE.0)call dinit(norder,zero,wksp1,1)
C
C     dealing with the source term:
C     this is distributed as using K_i^+ N^{-1} as the distribution matrix, i.e. LDA_2
C
C     1. Compute -N^{-1} = M^{-1}
C
      CALL DCOPY(MAXNOFVAR*MAXNOFVAR,SUM_K_NEG,1,SUM_K_POS,1) ! this is because SUM_K_NEG is later used to compute the approximate jacobian
      CALL DGETRI(NORDER,SUM_K_POS,MAXNOFVAR,IPIV,WKSP4,MAXNOFVAR,INFO) ! use WKSP4 as a workarray, overwriting its original content
!     WRITE(6,*)INFO
C
C     2. compute -N^{-1}*SOURCE and put into WR
C
      CALL DGEMV('No', NORDER, NORDER, MONE, SUM_K_POS, MAXNOFVAR,  ! -1 because SUM_K_POS is in fact -SUM_K_POS
     &SOURCE, 1, ZERO, WKSP4, 1)
C
C     distribute the convective term 
C
      DO 35 IVERT = 1,NOFVERT
          DO 33 IVAR = 1,NORDER
             NODRES(IVAR,IVERT) = PHIN(IVAR,IVERT) + ALFA(IVAR) *
     &       DPHI(IVAR,IVERT)
   33     CONTINUE
          IF(ICHECK.NE.0)
     +    call daxpy(norder,mone,nodres(1,ivert),1,wksp1,1)
   35 CONTINUE
caldo
#ifdef DEBUG
      CALL DINIT(NORDER,ZERO,DW,1)
#endif
caldo
C
C     add the source term: K^{+}_i N^{-1}*SOURCE
C
      DO IVERT = 1,NOFVERT
         CALL DGEMV('No', NORDER, NORDER, ONE, KPOS(1,1,IVERT),
     & MAXNOFVAR, WKSP4, 1, ONE, NODRES(1,IVERT), 1)
C
C    Here we check that what has been split among the nodes equals
C    the sum of the inviscid fluctuation and source term
C
#ifdef DEBUG
!
          CALL DAXPY(NORDER,ONE,NODRES(1,IVERT),1,DW,1) ! sum up the split fluctuations in DW
#endif
      ENDDO !  End Loop over nodes
#ifdef DEBUG
      CALL DAXPY(NORDER,ONE,CELRES,1,DW,1)
      CALL DAXPY(NORDER,ONE,SOURCE,1,DW,1)
      IF(DNRM2(NORDER,DW,1).GT.1.D-14)then
         write(6,*)'Smthg wrong DW @ ie =',ielem,(DW(I),i=1,norder)
      endif
#endif
C
C
C --------------- Debugging code starts here ---------------
C
      IF (ICHECK.NE.0) THEN
C
C     WKSP1 is the residual computed as \sum \Phi_j^{FCT}
C     it is compared with CELRES[1:NORDER] = - \sum K_j U_j
C
C     here the error seems far too large
C
          DO 12 IVAR = 1,NORDER
              TEMP = WKSP1(IVAR) - CELRES(IVAR)
              IF (DABS(TEMP).GT.1.D-15) THEN
                  WRITE (6,FMT=*) 'VKIB system scheme, elem ',IELEM,
     +              ' var # ',IVAR
                  WRITE (6,FMT=*)' computed ',WKSP1(IVAR),' "true" ',
     +              CELRES(IVAR),' E = ',TEMP
              ENDIF

   12     CONTINUE
      ENDIF
C
C --------------- Debugging code ends here ---------------
C
      IF (.NOT.MATRIX_ASSEMBLY) RETURN
C
C     Assembling the element stiffness matrix for the N scheme ..
C
      DO 31 J = 1,NOFVERT
C
C     Solve [\sum_K^{-}] Delta_j^{-} = K_j^{-}
C         Delta_j^{-} is overwritten onto K_j^{-}
C
          CALL DGETRS('N',NORDER,NORDER,SUM_K_NEG,MAXNOFVAR,IPIV,
     +                KNEG(1,1,J),MAXNOFVAR,INFO)
C
          DO 31 I = 1,NOFVERT
C
              IF (I.EQ.J) GOTO 31
C
C     C_{ij}^N = K_i^{+} Delta_j^{-} i neq j
C
              CALL DGEMM('N','N',NORDER,NORDER,NORDER,ONE,KPOS(1,1,I),
     +        MAXNOFVAR,KNEG(1,1,J),MAXNOFVAR,ZERO,STIFC(1,1,I,J),NDOF)
C
C     C_{ii} = - \sum_{j \neq i} C_{ij}
C                = - K_i^{+} [ \sum_{j \neq i} \Delta_j^- ] ...
C
              CALL MATSUB(STIFC(1,1,I,I),NDOF,STIFC(1,1,I,J),NDOF,
     +                    NORDER,NORDER)
C
   31 CONTINUE
C
C --------------- Debugging code starts here ---------------
C
      IF (ICHECK.NE.0) THEN
C
C        check that \sum_j \Delta_j^- = Identity matrix
C
          CALL DINIT(MAXNOFVAR**2,ZERO,WKSP3,1)
          DO 14 J = 1,NOFVERT
              CALL MATSUM(WKSP3,MAXNOFVAR,KNEG(1,1,J),
     +MAXNOFVAR,NORDER,NORDER)
   14     CONTINUE
          LFLAG = UNITMAT(WKSP3,NORDER,NORDER,MAXNOFVAR,1.D-14)
          IF (.NOT.LFLAG) THEN
              WRITE (6,FMT=*) IELEM,I
              CALL R8Mat_Print('General',' ',NORDER,NORDER,WKSP3,
     +             MAXNOFVAR,' Sum Delta_j^- equals identity ?',INFO)
caldo         PAUSE

          ENDIF
C
      include 'test2.inc'
C
      ENDIF
C
C --------------- Debugging code ends here ---------------
C
      RETURN

      END
