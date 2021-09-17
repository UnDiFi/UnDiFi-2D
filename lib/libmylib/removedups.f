      subroutine removedups(result,strings,n,nn)
      implicit none 
      integer n,nn
      integer result(*)
      integer strings(*)
      integer last_i,old_i
C
      result(1) = strings(1)

C   ' Copy the other items
!234567
      last_i = 1
      do old_i = 2,n 
        If(result(last_i) .NE. strings(old_i)) Then
            last_i = last_i + 1
            result(last_i) = strings(old_i)
        EndIf
      enddo
      nn = last_i
      return
      end

