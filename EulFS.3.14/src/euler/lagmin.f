       SUBROUTINE LAGMIN(GRADP,GRADD,GRADE,CHIC,SIG,KAPPAC,CHIAVG,KAVG)
C       
C      ##############################################################
C      #       Minimization problem for Lagrangian multiplier       #
C      #                                                            #
C      # input:                                                     # 
C      # GRADP = Pressure gradient                                  #
C      # GRADD = Gradient of chemical species density               #
C      # GRADE = Gradient of energy per unit of volume              #
C      # CHIC,KAPPAC = Chi and kappa approximation                  #
C      # SIG = Mean sigma (a^2)                                     #
C      #                                                            #
C      # output:                                                    # 
C      # CHIAVG,KAVG = Average of chi and kappa                     #
C      #                                                            #
C      ##############################################################
C
       INCLUDE 'paramt.h'
       INCLUDE 'plasma.h'

       INTEGER NMAXV,ISP,JSP,ID,JID,NID,NJD
       PARAMETER(NMAXV=NSP+1+NDIM)
       DOUBLE PRECISION GRADP(NDIM),GRADD(NSP,NDIM),GRADE(NDIM)
       DOUBLE PRECISION CHIC,KAPPAC,CHIAVG,KAVG,SIG,XIC(NSP),OMEGA
       DOUBLE PRECISION COEFM(NMAXV,NMAXV),KNOWT(NMAXV)
       DOUBLE PRECISION IPIVF,IPIVS,LDA,LDB,INFOF,INFOS
       CHARACTER*1 TRANS
C
C      Matrix of coefficients
C
       DO ISP = 1,NSP
          DO JSP = 1,NSP
             COEFM(ISP,JSP) = 0.d0
             IF(ISP.EQ.JSP)THEN
                COEFM(ISP,JSP) = 2/SIG**2
             ENDIF
           ENDDO 
           COEFM(ISP,IE) = 0.d0
           COEFM(ISP,IX) = GRADD(ISP,1)
           COEFM(ISP,IY) = GRADD(ISP,2)
           COEFM(IE,ISP) = 0.d0
           COEFM(IX,ISP) = GRADD(ISP,1)
           COEFM(IY,ISP) = GRADD(ISP,2)
           IF(NDIM.EQ.3)THEN
              COEFM(ISP,IZ) = GRAD(ISP,3)
              COEFM(IZ,ISP) = GRAD(ISP,3)
           ENDIF
       ENDDO        
C
       COEFM(IE,IE) = 2.d0
       DO ID = 1,NDIM
          NID = NSP + 1.d0 + ID
          COEFM(IE,NID) = GRADP(ID)
          COEFM(NID,IE) = GRADP(ID)
          DO JD = 1,NDIM
             NJD = NSP + 1.d0 + JD
             COEFM(NID,NJD) = 0.d0
          ENDDO 
       ENDDO
C
C      Vector of known terms
C
       OMEGAC = 1.d0/KAPPAC                   !omega = 1/k 
C
       DO ISP = 1,NSP
          XIC(ISP) = CHIC(ISP)/KAPPAC         !xi(isp) = chi(isp)/k
          KNOWT(ISP) = XIC(ISP)*2/SIG**2 
       ENDDO
C
       KNOWT(IE) = 2*OMEGAC
C
       DO ID = 1,NDIM
          NID = NSP + 1.d0 + ID 
          KNOWT(NID) = GRADE(ID)
       ENDDO
C
C      LU factorization
C
       LDA = NMAXV
       CALL DGETRF(NMAXV,NMAXV,COEFM,LDA,IPIVF,INFOF)
C
C      linear system solution 
C 
       TRANS = 'N'
       LDB = LDA
       CALL DGETRS(TRANS,NMAXV,NMAXV,COEFM,LDA,IPIVS,KNOWT,LDB,INFOS)
C
C      computing of averaged kappa and chi
C
       KAVG = 1/KNOWT(IE)
       DO ISP = 1,NSP
          CHIAVG(ISP) = KNOWT(ISP)*KAVG
       ENDDO     
     
       END
