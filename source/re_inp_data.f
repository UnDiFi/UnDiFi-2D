! Read the file input.dat containing the information concerning
! mesh generation in the proximity of the shocks shock integration
! additional hole points

      subroutine re_inp_data

      implicit none
      include 'paramt.h'

!     .. local scalars ..
      integer*4 i

!     open log file
      open(8,file='log/re_inp_data.log')

! open file input.dat containing the information concerning
!  1) mesh generation in the proximity of the shocks
!  2) shock integration
!  3) additional hole points
!  4) periodic boundarìes

      open(12,file='input.dat',status='old')
      write(8,*)' open file input.dat'
      read(12,*) eps        !< distance between the two shock faces
      read(12,*) sndmin     !< maximum nondimensional distance of the phantom nodes
      read(12,*) dxcell     !< lenght of shock edges
      read(12,*) shrelax    !< relax coefficent of shock point integration
      read(12,*) ibak       !< number of steps between the writing of a solution file and the next
      read(12,*) ga         !< specific heat ratio
      gm1=ga-1.0d+0

      read(12,*)naddholes   !< number of addition hole points
      do i=1,naddholes
       read(12,*)caddhole(1,i),caddhole(2,i) !< coordinates of addition hole points
      end do
      read(12,*)nprdbnd     !< number of periodic boundaries (allowed value 0=no periodic boundary 1 = one periodic boundary
      do i=1,nprdbnd
       read(12,*)prdbndclr(1,i), prdbndclr(2,i), prdbndclr(3,i) !< corresponding color pairs (indices 1 and 2)
                                                                !< of each periodic boundaries and equal coordinate
                                                                !< of periodic boundary (index 3 1=x 2=y)
      end do
      read(12,*)flt_dspeed  !< filter on discontinuity speeds [0,1.0] (0=disactive)
      read(12,*)imtf        !< iteration of mesh topology freezing (this value must be equal to one backup iteration)

      close(12)
      close(8)
      return
      end
