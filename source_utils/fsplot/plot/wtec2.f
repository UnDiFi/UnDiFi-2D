C
      subroutine WTecplotHeader(iunit,fheader,numq)
C
      IMPLICIT NONE
c  set max indicies in W and Y directions (may need to be increased)
      integer iunit,numq
c  initialize variables      
      character*(*) fheader
C
      CHARACTER FName*10,Text*14
      CHARACTER Variables*255
      INTEGER IERR
      INTEGER k,I,J
C234567890123456789012345678901234567890123456789012345678901234567890
      Variables(1:15) = 'VARIABLES = "Y"'
      k = 16
      do i = 1, numq
         write(Variables(k:k+8),FMT=220)i
         k = k+8
      enddo
  220 format(',"V(',I2,')"')
      write(6,FMT="(A)")Variables(1:k)
C
C     write header of Tecplot file
C
      WRITE(iunit,444)fheader,Variables(1:k)
  444 FORMAT('TITLE     = "',A30,'"',/,A)
C
      return
      end
      subroutine WriteXYZone(dfin,q,numq,npts)
C
      IMPLICIT NONE
c  set max indicies in W and Y directions (may need to be increased)
      integer numq
c  initialize variables      
      integer npts
      real*8 q(numq+1,*)
      character*(*) dfin
C
      CHARACTER FName*10,Text*14
      INTEGER nvt,IERR
      INTEGER k,I,J
C234567890123456789012345678901234567890123456789012345678901234567890
C
  446 FORMAT('ZONE T = "',A,'"',
     +/'I=',I3,', J=',I3,', K=1,F=POINT')
      WRITE(20,446)dfin,npts,1
         do 30 j=1,npts
            write (20,*)q(numq+1,j),(q(k,j),k=1,numq)
   30 continue
!     WRITE(6,333)'Closing',FName
      return
      end
C     ..
      subroutine WriteMacro(irake,xyz,ndim,npts)
      integer irake,npts
      double precision xyz(ndim,npts)
      double precision alpha
      CHARACTER MACRONAME*9
      integer kunit,i,j
      PARAMETER (alpha=305./71.7)
      DATA MACRONAME/"xh000.mcr"/
C
C     write a Tecplot macro for slicing through the mesh
C
      KUNIT = 22
      WRITE(MACRONAME(3:5),FMT="(I3.3)")IRAKE
      OPEN(UNIT=KUNIT,FILE=MACRONAME)
      WRITE(KUNIT,FMT=747)IRAKE,NPTS
  747 FORMAT('#!MC 800',/,'$!EXTRACTFROMPOLYLINE',/,
     +2X,'EXTRACTTHROUGHVOLUME=TRUE',/,
     +2X,'EXTRACTLINEPOINTSONLY=TRUE',/,
     +2X,'EXTRACTTOFILE=TRUE',/,2X,'FNAME="|BASEDIR|/',I3.3,'.dat"'/,
     +2X,'RAWDATA',/,I4)
C
      write(6,*)'Hello !' 
C
      DO 1 I = 1,NPTS
C
C ad hoc scaling
C
      xyz(1,i)= xyz(1,i)*alpha
      xyz(2,i)= xyz(2,i)*alpha
      xyz(3,i)=(xyz(3,i)*alpha)+3.
      WRITE(KUNIT,*)(xyz(j,i),j=1,ndim)
    1 CONTINUE
      CLOSE(UNIT=KUNIT)
      RETURN
      END
