NSYS_SCHEME                     Source Listing                  10-Aug-2002 18:29:04  DIGITAL Fortran 90 V5.1-594       Page 1
                                                                10-Aug-2002 18:02:03  NSys_scheme.f

	      1       SUBROUTINE NSYS_SCHEME(MATRIXSPLITTER,W,NODRES,DT,STIFC,NORDER,
	      2      +                       NOFVAR,NOFVERT,VCN,NDIM,DFGHDU,LDJ,CELRES,
	      3      +                       IELEM,TIMEIMPL)
	      4 C
	      5 C
	      6 C
	      7 C
	      8 C NORDER                IN Integer
	      9 C                       is the order of the system to solve for, i.e.
	     10 C                       the order of the matrix DFGHDU.
	     11 C IELEM                 IN Integer
	     12 C                       is the current element.
	     13 C DFGHDU(LDJ,*) is the Jacobian Matrix of the system.
	     14 C LDJ                   IN Integer
	     15 C                       is the leading dimension of DFGHDU.
	     16 C W(NORDER:NOFVERT)     stores by columns the NORDER variables of the
	     17 C                       NOFVERT vertices.
	     18 C LNODRES                   IN Integer
	     19 C                       is the leading dimension of Q.
	     20 C CELRES[1:2*NORDER]  OUT Real
	     21 C CELRES[1:NORDER]    stores the residual computed by the Matrix scheme
	     22 C                       as \sum K_j U_j (explicit part of the scheme)
	     23 C CELRES[NORDER+1:2*NORDER]
	     24 C                       stores the residual computed by the Matrix scheme
	     25 C                       as \sum C_{ij} U_j (implicit part of the scheme)
	     26 C MatrixSplitter        is the procedure used to compute the eigenvector
	     27 C                       decomposition of the matrix DFGHDU.
	     28 C DT                 is the nodal timestep.
	     29 C NODRES(NORDER,NOFVERT)    is the nodal residual.
	     30 C
	     31 C
	     32 C
	     33 C
	     34 C This routine computes the system N scheme on one tetrahedron
	     35 C
	     36       include 'paramt.h'
	     43       include 'flags.com'
	     44 C
	     45 C
	     46 C
	     47 C
	     48 C
	     49 C
	     50 C
	     51 C
	     52 C
	     53 C
	     54 C
	     55 C
	     56 C
	     57 C     .. Parameters ..
	     71       DOUBLE PRECISION ZERO,ONE
	     72       PARAMETER (ZERO=0.D0,ONE=1.D0)
	     73 C     ..
	     74 C     .. Scalar Arguments ..
	     75       INTEGER IELEM,LDJ,NDIM,NOFVERT,NORDER
	     76       LOGICAL TIMEIMPL

NSYS_SCHEME                     Source Listing                  10-Aug-2002 18:29:04  DIGITAL Fortran 90 V5.1-594       Page 2
                                                                10-Aug-2002 18:02:03  NSys_scheme.f

	     77 C     ..
	     78 C     .. Array Arguments ..
	     79       DOUBLE PRECISION DFGHDU(LDJ,*),NODRES(NOFVAR,NOFVERT),DT(NOFVERT),
	     80      +                 CELRES(*),STIFC(NOFVAR,NOFVAR,NOFVERT,NOFVERT),
	     81      +                 VCN(NDIM,NOFVERT),W(NOFVAR,NOFVERT)
	     82 C     ..
	     83 C     .. Subroutine Arguments ..
	     84       EXTERNAL MATRIXSPLITTER
	     85 C     ..
	     86 C     .. Local Scalars ..
	     87       DOUBLE PRECISION SUM,TEMP,TMP
	     88       INTEGER I,INFO,IROW,IVAR,IVERT,J,JCOL,NOFVAR
	     89       LOGICAL LFLAG
	     90 C     ..
	     91 C     .. Local Arrays ..
	     92       DOUBLE PRECISION KMAT(NMAX,NMAX),KNEG(NMAX,NMAX,VMAX),
	     93      +                 KPOS(NMAX,NMAX,VMAX),SUM_K_NEG(NMAX,NMAX),
	     94      +                 UNEG(NMAX),VLEFT(NMAX,NMAX),VRIGHT(NMAX,NMAX),
	     95      +                 WKSP1(NMAX),WKSP2(NMAX),WKSP3(NMAX,NMAX),
	     96      +                 WNEG(NMAX),WPOS(NMAX),WR(NMAX),DW(NMAX*VMAX)
	     97       INTEGER IPIV(NMAX)
	     98 C     ..
	     99 C     .. External Functions ..
	    100       DOUBLE PRECISION DNRM2
	    101       LOGICAL NULLMAT,UNITMAT
	    102       EXTERNAL DNRM2,NULLMAT,UNITMAT
	    103 C     ..
	    104 C     .. External Subroutines ..
	    105       EXTERNAL DAXPY,DGEMM,DGEMV,DGETRF,DGETRS,DINIT,DSCAL,MATSUB,
	    106      +         MATSUM,X04CAF
	    107 C     ..
	    108 C     .. Intrinsic Functions ..
	    109       INTRINSIC DABS
	    110 C     ..
	    111       DO 1 J = 1,NMAX
	    112           UNEG(J) = ZERO
	    113           DO 1 I = 1,NMAX
	    114               SUM_K_NEG(I,J) = ZERO
	    115     1 CONTINUE
	    116       DO 3 I = 1,2*NORDER
	    117           CELRES(I) = ZERO
	    118     3 CONTINUE
	    119 C
	    120 C --------------- Debugging code starts here ---------------
	    121 C
	    122       IF (ICHECK.NE.0) THEN
	    123           CALL DINIT(NMAX,ZERO,WKSP1,1)
	    124           CALL DINIT(NMAX,ZERO,WKSP2,1)
	    125           CALL DINIT(NMAX**2,ZERO,WKSP3,1)
	    126       ENDIF
	    127 C
	    128 C --------------- Debugging code ends here ---------------
	    129 C
	    130       DO 10 IVERT = 1,NOFVERT
	    131 C
	    132 C
	    133 C    The matrix is split into its positive and negative parts

NSYS_SCHEME                     Source Listing                  10-Aug-2002 18:29:04  DIGITAL Fortran 90 V5.1-594       Page 3
                                                                10-Aug-2002 18:02:03  NSys_scheme.f

	    134 C
	    135 C    Numerically or Analitically ..
	    136 C
	    137           CALL MATRIXSPLITTER(IELEM,NDIM,NORDER,VCN(1,IVERT),DFGHDU,LDJ,
	    138      +                        KMAT,KPOS(1,1,IVERT),KNEG(1,1,IVERT),
	    139      +                        VLEFT,VRIGHT,NMAX,WR,WPOS,WNEG,.TRUE.)
	    140 C
	    141 C    The negative jacobians are accumulated in SUM_K_NEG ..
	    142 C
	    143           DO 15 JCOL = 1,NORDER
	    144               DO 15 IROW = 1,NORDER
	    145                   SUM_K_NEG(IROW,JCOL) = SUM_K_NEG(IROW,JCOL) +
	    146      +                                   KNEG(IROW,JCOL,IVERT)
	    147    15     CONTINUE
	    148 C
	    149 C    Computes KNEG(IVERT)*U(IVERT) and adds it to UNEG ..
	    150 C
	    151 cblas    CALL DGEMV('N',NORDER,NORDER,ONE,KNEG(1,1,IVERT),NMAX,W
	    152 cblas+   (1,IVERT),1,ONE,UNEG,1)
	    153 C
	    154           DO 11 I = 1,NORDER
	    155               SUM = ZERO
	    156               DO 13 J = 1,NORDER
	    157                   SUM = SUM + KNEG(I,J,IVERT)*W(J,IVERT)
	    158    13         CONTINUE
	    159               UNEG(I) = UNEG(I) + SUM
	    160    11     CONTINUE
	    161 C
	    162 C
	    163 C --------------- Debugging code starts here ---------------
	    164 C
	    165           IF (ICHECK.NE.0) THEN
	    166 C
	    167 C    the residual is computes as Sum_j K(j) * U(j) ..
	    168 C
	    169               CALL DGEMV('N',NORDER,NORDER,ONE,KMAT,NMAX,W(1,IVERT),1,
	    170      +                   ONE,CELRES,1)
	    171           ENDIF
	    172 C
	    173 C --------------- Debugging code ends here ---------------
	    174 C
	    175 C    Timestep ..
	    176 C
	    177           DO 34 I = 1,NORDER
	    178               DT(IVERT) = DT(IVERT) + WPOS(I)
	    179    34     CONTINUE
	    180 C
	    181    10 CONTINUE
	    182 C
	    183 C    Finds the generalized inflow point
	    184 C       solving [\sum_j K_j^-] U_{-} = \sum_j K_j^- U_j
	    185 C
	    186 C    LU factorization ..
	    187 C
	    188 cnag  CALL F07ADF(NORDER,NORDER,SUM_K_NEG,NMAX,IPIV,INFO)
	    189       CALL DGETRF(NORDER,NORDER,SUM_K_NEG,NMAX,IPIV,INFO)
	    190 C

NSYS_SCHEME                     Source Listing                  10-Aug-2002 18:29:04  DIGITAL Fortran 90 V5.1-594       Page 4
                                                                10-Aug-2002 18:02:03  NSys_scheme.f

	    191       IF (INFO.GT.0) THEN
	    192           WRITE (6,FMT=99999) IELEM
	    193 
	    194 99999     FORMAT (5X,'Matrix SUM_K_NEG is singular in IELEM = ',I6)
	    195 *
	    196           DO 9 IVERT = 1,NOFVERT
	    197               WRITE (6,FMT=*) '    Vertice # ',IVERT
	    198               CALL X04CAF('G',' ',NORDER,NORDER,KPOS(1,1,IVERT),NMAX,
	    199      +                    'K(+) ',INFO)
	    200               CALL X04CAF('G',' ',NORDER,NORDER,KNEG(1,1,IVERT),NMAX,
	    201      +                    'K(-) ',INFO)
	    202     9     CONTINUE
	    203           CALL X04CAF('G',' ',NORDER,NORDER,SUM_K_NEG,NMAX,
	    204      +                'SUM_j K(-) ',INFO)
	    205           STOP
	    206 
	    207       ENDIF
	    208 C
	    209 C    solution ..
	    210 C
	    211 cnag  CALL F07AEF('N',NORDER,1,SUM_K_NEG,NMAX,IPIV,UNEG,NMAX,INFO)
	    212       CALL DGETRS('N',NORDER,1,SUM_K_NEG,NMAX,IPIV,UNEG,NMAX,INFO)
	    213 C
	    214 C Loops over nodes
	    215 C
	    216       DO 30 IVERT = 1,NOFVERT
	    217           I = (IVERT-1)*NORDER
	    218 !         DO 35 IVAR = 1,NORDER
	    219 !             DW(I+IVAR) = W(IVAR,IVERT) - UNEG(IVAR)
	    220 !  35     CONTINUE
	    221 C
	    222 C    NODRES[INODE] = - KPOS(IVERT)*[U(IVERT)-UNEG]
	    223 C
	    224 cblas    CALL DGEMV('N',NORDER,NORDER,ONE,KPOS(1,1,IVERT),NMAX, W
	    225 cblas+   (1,IVERT),1,ZERO,WKSP,1)
	    226 C
	    227           DO 23 IVAR = 1,NORDER
	    228               SUM = ZERO
	    229               TMP = ZERO
	    230               DO 21 J = 1,NORDER
	    231                   TEMP = KPOS(IVAR,J,IVERT)
	    232                   SUM = SUM + TEMP * (W(J,IVERT)-UNEG(J))
	    233                   TMP = TMP + TEMP * W(J,IVERT)
	    234    21         CONTINUE
	    235               NODRES(IVAR,IVERT) = -SUM
	    236    23     CONTINUE
	    237 C
	    238 C --------------- Debugging code starts here ---------------
	    239 C
	    240           IF (ICHECK.NE.0) THEN
	    241 C
	    242 C    The residual is computed as Sum_j KPOS(j) * [U(j)-Uin]
	    243 C       and stored in WKSP1
	    244 C
	    245 
	    246               CALL DAXPY(NORDER,-ONE,NODRES(1,IVERT),1,WKSP1,1)
	    247           ENDIF

NSYS_SCHEME                     Source Listing                  10-Aug-2002 18:29:04  DIGITAL Fortran 90 V5.1-594       Page 5
                                                                10-Aug-2002 18:02:03  NSys_scheme.f

	    248 C
	    249 C --------------- Debugging code ends here ---------------
	    250 C
	    251    30 CONTINUE
	    252 C
	    253 C  End Loops over nodes
	    254 C
	    255 C
	    256 C --------------- Debugging code starts here ---------------
	    257 C
	    258       IF (ICHECK.NE.0) THEN
	    259 C
	    260 C     WKSP1 is the residual computed as \sum \Phi_j^{N}
	    261 C     it is compared with CELRES[1:NORDER] = - \sum K_j U_j
	    262 C
	    263           DO 12 IVAR = 1,NORDER
	    264               TEMP = WKSP1(IVAR) - CELRES(IVAR)
	    265               IF (DABS(TEMP).GT.1.D-15) THEN
	    266                   WRITE (6,FMT=*) 'N system scheme, elem ',IELEM,
	    267      +              ' var # ',IVAR,' computed ',WKSP1(IVAR),' "true" ',
	    268      +              CELRES(IVAR)
	    269               ENDIF
	    270 
	    271    12     CONTINUE
	    272       ENDIF
	    273 C
	    274 C --------------- Debugging code ends here ---------------
	    275 C
	    276       IF (.NOT.TIMEIMPL) RETURN
	    277 C
	    278 C     Assembling the element stiffness matrix for the N scheme ..
	    279 C
	    280       DO 31 J = 1,NOFVERT
	    281 C
	    282 C     Solve [\sum_K^{-}] Delta_j^{-} = K_j^{-}
	    283 C         Delta_j^{-} is overwritten onto K_j^{-}
	    284 C
	    285           CALL DGETRS('N',NORDER,NORDER,SUM_K_NEG,NMAX,IPIV,KNEG(1,1,J),
	    286      +                NMAX,INFO)
	    287 C
	    288           DO 31 I = 1,NOFVERT
	    289 C
	    290               IF (I.EQ.J) GOTO 31
	    291 C
	    292 C     C_{ij}^N = K_i^{+} Delta_j^{-} i neq j
	    293 C
	    294               CALL DGEMM('N','N',NORDER,NORDER,NORDER,ONE,KPOS(1,1,I),
	    295      +             NMAX,KNEG(1,1,J),NMAX,ZERO,STIFC(1,1,I,J),NOFVAR)
	    296 C
	    297 C     C_{ii} = - \sum_{j \neq i} C_{ij}
	    298 C                = - K_i^{+} [ \sum_{j \neq i} \Delta_j^- ] ...
	    299 C
	    300               CALL MATSUB(STIFC(1,1,I,I),NOFVAR,STIFC(1,1,I,J),NOFVAR,
	    301      +                    NORDER,NORDER)
	    302 C
	    303    31 CONTINUE
	    304 C

NSYS_SCHEME                     Source Listing                  10-Aug-2002 18:29:04  DIGITAL Fortran 90 V5.1-594       Page 6
                                                                10-Aug-2002 18:02:03  NSys_scheme.f

	    305 C --------------- Debugging code starts here ---------------
	    306 C
	    307       IF (ICHECK.NE.0) THEN
	    308 C
	    309 C        check that \sum_j \Delta_j^- = Identity matrix
	    310 C
	    311           CALL DINIT(NMAX**2,ZERO,WKSP3,1)
	    312           DO 14 J = 1,NOFVERT
	    313               CALL MATSUM(WKSP3,NMAX,KNEG(1,1,J),NMAX,NORDER,NORDER)
	    314    14     CONTINUE
	    315           LFLAG = UNITMAT(WKSP3,NORDER,NORDER,NMAX,1.D-14)
	    316           IF (.NOT.LFLAG) THEN
	    317               WRITE (6,FMT=*) IELEM,I
	    318               CALL X04CAF('General',' ',NORDER,NORDER,WKSP3,NMAX,
	    319      +                    ' Sum Delta_j^- equals identity ?',INFO)
	    320               PAUSE
	    321 
	    322           ENDIF
	    323 C
	    324 C     Check that \sum_{j} C_{i,j} = 0
	    325 C
	    326           DO 16 I = 1,NOFVERT
	    327               CALL DINIT(NMAX*NMAX,ZERO,WKSP3,1)
	    328               DO 17 J = 1,NOFVERT
	    329                   CALL MATSUM(WKSP3,NMAX,STIFC(1,1,I,J),NOFVAR,NORDER,
	    330      +                        NORDER)
	    331    17         CONTINUE
	    332               LFLAG = NULLMAT(WKSP3,NORDER,NORDER,NMAX,1.D-14)
	    333               IF (.NOT.LFLAG) THEN
	    334                   WRITE (6,FMT=*) IELEM,I
	    335                   CALL X04CAF('General',' ',NORDER,NORDER,WKSP3,NMAX,
	    336      +                        ' Sum_j C_{i,j} equals the null matrix ?',
	    337      +                        INFO)
	    338                   PAUSE
	    339 
	    340               ENDIF
	    341 
	    342    16     CONTINUE
	    343 C
	    344           GOTO 69
	    345 C
	    346 C        print the element stiffness matrix
	    347 C
	    348           DO 19 J = 1,NOFVERT
	    349               DO 19 I = 1,NOFVERT
	    350                   WRITE (6,FMT=*) 'IELEM = ',IELEM,'I,J = ',I,J
	    351                   CALL X04CAF('General',' ',NORDER,NORDER,
	    352      +                        STIFC(1,1,I,J),NOFVAR,
	    353      +                        ' Element stiffness matrix ',INFO)
	    354    19     CONTINUE
	    355           PAUSE
	    356 
	    357    69     CONTINUE
	    358 C
	    359 C     CELRES[NORDER+1:2*NORDER] stores the residual 
	    360 C     in the "V" variable computed as - \sum_{i,j} C_{i,j} U_j
	    361 C

NSYS_SCHEME                     Source Listing                  10-Aug-2002 18:29:04  DIGITAL Fortran 90 V5.1-594       Page 7
                                                                10-Aug-2002 18:02:03  NSys_scheme.f

	    362           DO 22 J = 1,NOFVERT
	    363               DO 22 I = 1,NOFVERT
	    364                   CALL DGEMV('N',NORDER,NORDER,-ONE,STIFC(1,1,I,J),
	    365      +                       NOFVAR,W(1,J),1,ONE,CELRES(NORDER+1),1)
	    366    22     CONTINUE
	    367 C
	    368 C
	    369 C
	    370           LFLAG = .FALSE.
	    371           DO 27 IVAR = 1,NORDER
	    372               TEMP = CELRES(IVAR) - CELRES(NORDER+IVAR)
	    373               IF (DABS(TEMP).GT.1.D-14) LFLAG = .TRUE.
	    374    27     CONTINUE
	    375           IF (LFLAG .EQV. .TRUE.) THEN
	    376               WRITE (6,FMT=*) 'N system scheme, elem ',IELEM
	    377               WRITE (6,FMT=*) 'var #     implicit        "true"'
	    378               DO 25 IVAR = 1,NORDER
	    379                   WRITE (6,FMT=*) IVAR,CELRES(NORDER+IVAR),
	    380      +              CELRES(IVAR)
	    381    25         CONTINUE
	    382           ENDIF
	    383 C
	    384       ENDIF
	    385 C
	    386 C --------------- Debugging code ends here ---------------
	    387 C
	    388       RETURN
	    389 
	    390       END

NSYS_SCHEME                     Source Listing                  10-Aug-2002 18:29:04  DIGITAL Fortran 90 V5.1-594       Page 8
                                                                10-Aug-2002 18:02:03  NSys_scheme.f



PROGRAM SECTIONS
    Name				 Bytes   Attributes
  1 $BSS$                                 4320 NOPIC CON REL LCL NOSHR NOEXE   RD   WRT OCTA
  2 $DATA$                                  76 NOPIC CON REL LCL NOSHR NOEXE   RD   WRT OCTA
  3 $CODE$                                7476   PIC CON REL LCL   SHR   EXE NORD NOWRT OCTA
  4 cflags_                                 32 NOPIC OVR REL GBL NOSHR NOEXE   RD   WRT OCTA
  5 .lita                                  176   PIC CON REL LCL NOSHR NOEXE   RD NOWRT QUAD
  6 .lit4                                   44   PIC CON REL LCL NOSHR NOEXE   RD NOWRT QUAD
  7 .lit8                                   64   PIC CON REL LCL NOSHR NOEXE   RD NOWRT QUAD
  8 .literal                               233   PIC CON REL LCL NOSHR NOEXE   RD NOWRT QUAD
  9 .xdata                                   8 NOPIC CON REL LCL NOSHR NOEXE NORD NOWRT QUAD
 10 .pdata                                   8 NOPIC CON REL LCL NOSHR NOEXE NORD NOWRT QUAD

    Total Space Allocated                12437


ENTRY POINTS

 Address      Name        
                          
  3-00000000  nsys_scheme_


VARIABLES
 Address      Type  Name           Address      Type  Name           Address      Type  Name         
                                                                                                     
  4-00000010  I*4   DECOMP          1-00000030  I*4   J                   **      I*4   NORDER       
  1-00000044  I*4   I               1-0000002C  I*4   JCOL            4-0000001C  L*4   PERIODIC_MESH
  4-00000008  I*4   ICHECK          4-00000004  I*4   JSCHEME         1-00000010  R*8   SUM          
      **      I*4   IELEM           4-0000000C  I*4   KAN             1-00000008  R*8   TEMP         
  1-00000040  I*4   INFO                **      I*4   LDJ                 **      L*4   TIMEIMPL     
  1-0000003C  I*4   IROW            1-00000028  L*4   LFLAG           1-00000000  R*8   TMP          
  4-00000000  I*4   ISCHEME             **      I*4   NDIM            4-00000018  L*4   TURBULENT    
  1-00000038  I*4   IVAR                **      I*4   NOFVAR          4-00000014  I*4   WBCTYP       
  1-00000034  I*4   IVERT               **      I*4   NOFVERT                                        


ARRAYS

 Address      Type  Name             Bytes  Dimensions

      **      R*8   CELRES              **  (*)
      **      R*8   DFGHDU              **  (*, *)
      **      R*8   DT                  **  (*)
  1-00000A80  R*8   DW                 192  (24)
  1-00000168  I*4   IPIV                24  (6)
  1-00000FC0  R*8   KMAT               288  (6, 6)
  1-00000600  R*8   KNEG              1152  (6, 6, 4)
  1-00000180  R*8   KPOS              1152  (6, 6, 4)
      **      R*8   NODRES              **  (*, *)
      **      R*8   STIFC               **  (*, *, *, *)
  1-00000EA0  R*8   SUM_K_NEG          288  (6, 6)
  1-00000138  R*8   UNEG                48  (6)
      **      R*8   VCN                 **  (*, *)
  1-00000D80  R*8   VLEFT              288  (6, 6)

NSYS_SCHEME                     Source Listing                  10-Aug-2002 18:29:04  DIGITAL Fortran 90 V5.1-594       Page 9
                                                                10-Aug-2002 18:02:03  NSys_scheme.f

  1-00000C60  R*8   VRIGHT             288  (6, 6)
      **      R*8   W                   **  (*, *)
  1-00000108  R*8   WKSP1               48  (6)
  1-000000D8  R*8   WKSP2               48  (6)
  1-00000B40  R*8   WKSP3              288  (6, 6)
  1-000000A8  R*8   WNEG                48  (6)
  1-00000078  R*8   WPOS                48  (6)
  1-00000048  R*8   WR                  48  (6)

LABELS
 Address      Label  Address      Label  Address      Label  Address      Label  Address      Label  Address      Label
                                                                                                                       
      **      1           **      11          **      15          **      21      3-00001AF4  27          **      69   
      **      3       3-00001058  12      3-00001874  16          **      22      3-00000DF4  30                       
      **      9           **      13          **      17          **      23      3-000013B8  31                       
  3-00000804  10          **      14          **      19          **      25          **      34                       


FUNCTIONS AND SUBROUTINES REFERENCED
  Type  Name                
                            
         MATRIXSPLITTER     

NSYS_SCHEME                     Source Listing                  10-Aug-2002 18:29:04  DIGITAL Fortran 90 V5.1-594       Page 10
                                                                10-Aug-2002 18:02:03  NSys_scheme.f

COMPILER OPTIONS BEING USED

  no   -align commons                   no   -align dcommons
       -align records                   no   -align rec1byte
  no   -align rec2byte                  no   -align rec4byte
  no   -align rec8byte                       -altparam
       -arch generic                         -assume accuracy_sensitive
  no   -assume bigarrays                no   -assume byterecl
  no   -assume dummy_aliases            no   -assume minus0
       -assume underscore                    -assume source_include
       -assume zsize                    no   -automatic
       -call_shared                     no   -check bounds
  no   -check format                    no   -check omp_bindings
  no   -check output_conversion         no   -check overflow
       -check power                     no   -check underflow
       -convert native                       -double_size 64
  no   -d_lines                              -error_limit 30
  no   -extend_source                   no   -f66
       -fixed                           no   -fpconstant
       -fpe0                                 -fprm nearest
       -g2                                   -granularity quadword
  no   -hpf_matmul                      no   -intconstant
       -integer_size 32                 no   -ladebug
  no   -machine_code                         -math_library accurate
  no   -module                               -names lowercase
       -nearest_neighbor                no   -nowsf_main
  no   -non_shared                      no   -noinclude
       -numnodes 0                           -O0
       -inline none                     no   -transform_loops
  no   -pipeline                             -speculate none
       -tune generic                         -unroll 0
  no   -pad_source                           -parallel manual
  no   -pg                                   -real_size 32
  no   -recursive                            -reentrancy none
       -shadow_width 0                  no   -shared
  no   -show include                         -show map
  no   -show wsfinfo                    no   -show hpf_all
  no   -show hpf_punt                   no   -show hpf_nearest
  no   -show hpf_comm                   no   -show hpf_temps
  no   -show hpf_indep                  no   -show hpf_dev
  no   -show hpf_default                no   -std
  no   -synchronous_exceptions          no   -syntax_only
  no   -vms                                  -warn alignments
  no   -warn argument_checking          no   -warn declarations
       -warn general                         -warn granularity
  no   -warn truncated_source                -warn uncalled
       -warn uninitialized                   -warn usage
       -warning_severity warning        no   -wsf
  no   -fuse_xref

       -I    path : ../../include/,/usr/lib/cmplrs/hpfrtl/,/usr/include/
       -V    filename : NSys_scheme.l
       -o    filename : NSys_scheme.o

COMPILER: DIGITAL Fortran 90 V5.1-594-3882K

NSYS_SCHEME                     Source Listing                  10-Aug-2002 18:29:04  DIGITAL Fortran 90 V5.1-594       Page 11
                                                                10-Aug-2002 18:02:03  NSys_scheme.f

COMPILATION STATISTICS

  CPU time:          0.00 seconds
  Elapsed time:      0.00 seconds
  Pagefaults:           0
  I/O Count:            0
