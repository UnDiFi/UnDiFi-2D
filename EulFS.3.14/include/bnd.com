      INTEGER ICOLOR(0:NCOLOR,3),MCOLOR(0:NCOLOR),IBGN(MBODIES),
     +        IEND(MBODIES),IMUNIT(0:NCOLOR),IFUNIT(0:NCOLOR),
     +        IBFLX(NCOLOR),NBFLX(2)
      DOUBLE PRECISION SCOLOR(0:NCOLOR),VISCF(3,0:MBODIES),
     +PRESF(3,0:MBODIES),CFLUX(0:NCOLOR),VOLTAGE(0:NCOLOR,-1:1)
      CHARACTER*20 CBTYPE(0:NCOLOR)
C
C     ICOLOR is the boundary type corresponding to a given color
C     i.e.   0 <= ICOLOR(i) <= NBTYPE
C     ICOLOR(*,1) is for the mean flow equations
C     ICOLOR(*,2) is for the turbulence transport equations
C     ICOLOR(*,3) is for Laplace's (Poisson's) equations
C     MCOLOR is the # of boundary faces for a given color
C     SCOLOR is the total surface of the boundary of a given color
C     IMUNIT is the unit number used to open files associated with each color
C     IBFLX(NBFLX(1)) is a list of boundary segments where a flux should be computed
C     NBFLX(1) is the nof boundary patches where the flux is prescribed
C     NBFLX(2) is the nof boundary edges where the flux is prescribed
C     VOLTAGE(0:NCOLOR,-1:1) is the voltage (non-dimensional) that is applied
C              on the boundary patches, the second index refers to the time-level
C              0 = current, 1 = next, -1 = old
C
      COMMON /bnd_i4/ ICOLOR,MCOLOR,
     &                IBGN,IEND,IMUNIT,IFUNIT,IBFLX,NBFLX
      COMMON /bnd_r8/ SCOLOR,VISCF,PRESF,CFLUX,VOLTAGE
      COMMON /bnd_ch/ CBTYPE
