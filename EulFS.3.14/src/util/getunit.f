      subroutine get_unit ( iunit )
C
C*******************************************************************************
C
C! GET_UNIT returns a free FORTRAN unit number.
C
C
C  Discussion:
C
C    A "free" FORTRAN unit number is an integer between 1 and 99 which
C    is not currently associated with an I/O device.  A free FORTRAN unit
C    number is needed in order to open a file with the OPEN command.
C
C  Modified:
C
C    02 March 1999
C
C  Author:
C
C    John Burkardt
C
C  Parameters:
C
C    Output, integer IUNIT.
C
C    If IUNIT = 0, then no free FORTRAN unit could be found, although
C    all 99 units were checked (except for units 5 and 6).
C
C    Otherwise, IUNIT is an integer between 1 and 99, representing a
C    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
C    are special, and will never return those values.
C
      implicit none
C
C234567
      integer i
      integer ios
      integer iunit
      logical lopen

      iunit = 0

      do i = 10, 99
c
          if ( i .NE. 5 .and. i .NE. 6 .AND. i .NE. 9) then

          inquire ( unit = i, opened = lopen, iostat = ios )

             if ( ios .EQ. 0 ) then
                if ( .not. lopen ) then
                  iunit = i
                  return
                end if
             end if

          end if

      end do

      return
      end
