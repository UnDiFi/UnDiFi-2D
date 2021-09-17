      subroutine subyy(icn,nofvert,klist,vlist,lda,nlist,pcn)
      implicit none
c
c     $Id: subyy.f,v 1.3 2005/09/22 08:17:02 abonfi Exp $     
c
c     locates a meshpoint among those in the list of
c     inflow nodes
c     in the current release this also includes
c     periodic nodes
c
      integer nofvert,nlist,lda
      integer icn(nofvert),klist(nlist)
      double precision vlist(lda,nlist),pcn(lda,nofvert)
      integer ipos,last,ivert,nerr,iopt,k
      character*55 ERRMSG
c
      do 1 ivert= 1, nofvert-1
c
         call binsrc(icn(ivert),klist,nlist,ipos,last)
         if(ipos.EQ.0)THEN
            WRITE(ERRMSG(1:55),FMT=300)icn(ivert)
            WRITE(6,FMT=*)(icn(k),k=1,nofvert)
            NERR = 4
            IOPT = 1
            CALL SETERR(ERRMSG,56,NERR,IOPT)
         ENDIF
         do 1 k = 1,lda
            pcn(k,ivert) = vlist(k,ipos)
    1 continue
      return
  300 FORMAT('SUBYY CANNOT LOCATE NODE ',I6.6,
     +       ' AMONG THOSE IN THE LIST')
      end
