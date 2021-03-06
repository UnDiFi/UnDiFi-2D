      INTEGER     LCORG,LZROE,LTURB,LCELNOD,LCELFAC,LCELCEL, ! 6
     1 LFACNOR,LBNDFAC,LVOL,LNODCOD,LSKINF,LHEAT,LTD,LTZX,LTTD,LCHPSI, ! 16
     2 LPMAP,LFREE,LMEDIAN,LCLRC,LCLA,LCLJA,LCLIA,LCLDEG,LCLZB,
     2 LPROBE(2),LBNDFLX(2),LXYZDOT ! 30
      COMMON/NLOC/LCORG,LZROE,LTURB,LCELNOD,LCELFAC,LCELCEL,
     1 LFACNOR,LBNDFAC,LVOL,LNODCOD,LSKINF,LHEAT,LTD,LTZX,LTTD,LCHPSI,
     2 LPMAP,LFREE,LMEDIAN,LCLRC,LCLA,LCLJA,LCLIA,LCLDEG,LCLZB,LPROBE,
     3 LBNDFLX,LXYZDOT
C
C     Common block for POINTERS 
C
C     LCORG   ---> CORG(1:NDIM,1:NPOIN)
C     LZROE   ---> ZROE(1:NOFVAR,1:NPOIN)
C     LCELNOD ---> ICELNOD(1:NOFVERT,1:NELEM)
C     LCELFAC ---> ICELFAC(1:NOFVERT,1:NELEM)
C     LCELCEL ---> ICELCEL(1:NOFVERT,1:NELEM)
C     LFACNOR ---> FACNOR(1:NDIM,1:NFACE)
C     LBNDFAC ---> IBNDFAC(1:3,1:NBFACE)
C     LVOL    ---> VOL(1:NELEM)
C     LNODCOD ---> INODCOD(1:NPOIN)
C     LSKINF  ---> SKIN(1:*) Skin friction coefficient on wall edges
C     LHEAT   ---> HEAT(1:*) Heat flux coefficient on wall edges
C     LTD     ---> TD(1:NPOIN) Nearest wall distance
C     LTZX    ---> unused
C     LPMAP   ---> mapping for periodic nodes
C     LPTOT   ---> map of inlet total pressure
C     LCLRC   ---> RCLINE(1:NCL) radial position of c-lines
C     LCLA    ---> CLA(4,1:NNR) area coordinates and theta position of c-lines
C     LCLJA   ---> CLJA(3,1:NNR) 3 meshpoints surrounding the j-th node of a c-lines
C     LCLIA   ---> CLIA(1:NCL+1) CSR style linked list for the c-lines
C     LCLZB   ---> CLZB(NOFVAR,1:NCL) CSR style linked list of c-lines
C     LCLDEG  ---> CLDEG(1:NCL) # entries for each c-line = IA(I+1)-IA(I)
C                  required in the parallel case
C     LPROBE(1)--> IELEM(1:NPROBE) element to which the probed (x,y,z) belongs to
C     LPROBE(2)--> WEIGHTS(1:NOFVERT,1:NPROBE) weigths for probe
C     LBNDFLX(1) --> BFLX(1:NOFVAR,*) NOFVAR components of the prescribed flux
C     LBNDFLX(2) --> IBFLX(*) addresses the entry of IBNDPTR(1:3,1:NBFAC)
C     LXYZDOT    --> x,y,z components of the nodal grid velocity
C
C
C	********** BE CAREFUL **********
C
C	If you modify the NLOC common, you have to modify the
C	same common in the routine "blockdata" where it is NOT included. 
C
