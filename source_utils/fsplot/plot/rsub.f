c
      subroutine rsub(npdes,nnode,q)

      implicit none
      integer*4 npdes,nnode,ipde,in,i,j
      real*8    q(npdes,nnode)
      real*8 ux,uy,uz,pres,kine,h,z1,zk
      parameter(zk=1.4)

      open (unit=3,file='hydra2aldo',form='unformatted',
     &      status='old')

      do in=1,nnode
        do ipde=1,npdes
          read(3,*) q(ipde,in) 
        end do
      end do

c     All variables are nondimensionalized.
c     q(1,in) density
c     q(2,3,4,in) 3 velocity components in Cartesian coordinates
c     (I think q(4,in) could be safely set to 0).
c     q(5,in) static pressure
c     q(6,in) turbulence-related variable. Set it to 0.001

c     Nondimensionalization

c     reference length Lref = 1 meter
c     reference press = 1.013e5
c     reference density (rref)  = 1.226
c     reference velocity uref = sqrt(pref/rref)
c     reference time = Lref / uref
c
c

      do in=1,nnode
         z1 = sqrt(q(1,in))
         ux = q(2,in)
         uy = q(3,in)
         uz = q(4,in)
         pres = q(5,in)
         kine = 0.5*( ux*ux+uy*uy+uz*uz )
         h = zk/(zk-1.)*q(5,in)/q(1,in) + kine
         q(3,in) = z1 * ux
         q(4,in) = z1 * uy
         q(5,in) = z1 * uz
         q(2,in) = z1 * h
      end do

      close(unit=3)

      return
      end
