module mod_setbndrynodeptr
! setbndrynodeptr's nodptr/ia/ja/iclr dummies are allocatable,
! intent(out) -- Fortran requires an explicit interface to call such a
! procedure correctly (the compiler must pass a full array descriptor,
! not just a base address, so the callee can allocate the actual). Its
! sole caller, readmesh.f90, previously declared it via a plain
! `external setbndrynodeptr` (implicit interface) instead -- a real
! standard violation that gfortran's -fcheck=all catches at runtime as
! "Allocatable actual argument 'mesh' is not allocated" at the call
! site (confirmed present since Phase 3.1's first commit, unrelated to
! any Phase 3.2 special-point work: readmesh.f90/setbndrynodeptr.f90
! are both Phase 2 files). This module supplies that explicit
! interface; setbndrynodeptr's own body is moved here verbatim from
! setbndrynodeptr.f90 (which still holds myroutine/findcolours/
! setbndrynodelist/nearby/check -- none of those have allocatable
! dummies, so calling them via implicit interface remains safe and is
! left untouched).
  implicit none(type, external)
  private
  public :: setbndrynodeptr

contains

!> @bndfac[in] bndry face pointer read from the poly file, bndfac(3,nbfac)
!> @nodcod[in] nodal pointer read from the node or poly file, nodcod(npoin)
!> @nbfac[in] nof bndry faces
!> @npoin[in] nof vertices
!> @nbpoin[out] nof boundary vertices
!> @nodptr[out] pointer for boundary vertices, nodptr(nbpoin,3)
!> @ia[out] ia(1:nclr+1)
!> @ja[out] ja(1:nnzr)
!> @iclr[out] iclr(1:nclr) colour of the nclr patches
  subroutine setbndrynodeptr(bndfac, nodcod, nbfac, npoin,&
  &nbpoin, nodptr, ia, ja, iclr, nclr)
!
    use mod_kinds, only: wp, i4
!
!     ndim  is the space dimension =2
!     nvt = ndim+1 is the number of vertices
!
!     npshmax     : max number of shock points.
!     neshmax     : max number of shock elements
!
!     .. scalar arguments ..
    integer(i4), intent(in) :: nbfac, npoin
    integer(i4), intent(out) :: nbpoin, nclr

!     .. array arguments ..
    integer(i4), intent(in) :: bndfac(3, nbfac), nodcod(npoin)
    integer(i4), allocatable, intent(out) :: nodptr(:, :), ia(:), ja(:), iclr(:)

!     .. local scalars ..
    integer(i4) ipoin, k, nnzr

!     .. local arrays ..
    integer(i4), parameter :: maxpatches = 50
    integer(i4) ic(0:maxpatches) ! number of bndry gridpoints within patch coloured i, 0 <=i<= maxpatches
    logical closed(0:maxpatches) ! .true. if the boundary patch is closed (i.e. a profile)
    integer(i4), allocatable :: iwork(:)
!     ..
!     .. external subroutines ..
    external findcolours, myroutine, setbndrynodelist
!
!     count the nof boundary vertices
!
    nbpoin = 0
    do ipoin = 1, npoin
      if (nodcod(ipoin) .gt. 0) nbpoin = nbpoin + 1
    end do
    write (6, *) 'setbndrynodeptr: found ', nbpoin, ' boundary points'
    allocate (nodptr(nbpoin, 3), source=0_i4)
    allocate (iwork(2*nbpoin), source=0_i4)
    call myroutine(nodcod, bndfac, nbfac, nodptr, npoin, nbpoin, iwork)
    deallocate (iwork) ! workspace
    call findcolours(bndfac, nodptr, ic, closed, maxpatches, nbfac, nbpoin, nclr)
    write (6, *) 'setbndrynodeptr: found ', nclr, ' colours or bndry patches'
!
    allocate (ia(nclr + 1), source=0_i4) ! pointer in csr format
    allocate (iclr(nclr), source=0_i4) ! colour of the k-th boundary patch
!
!     we use two pointers ia(1:nclr+1) and ja(1:nnzr) to address the bndry gridpoints
!     gridpoints belonging to boundary (or patch) i, where 1<=i<=nclr
!     are stored in ja(j), j=jbgn,jend where
!     jbgn = ia(i), jend = ia(i+1)-1
!     the total nof entries in ja is nnzr = ia(nclr+1)-ia(1);
!     note that ia(1) = 1
!
!     therefore, in order to address all gridpoints of patch coloured 3, we do:
!
!     jbgn = ia(3)
!     jend = ia(4) -1
!     do j = jbgn,jend
!        ipoin = ja(j) ! global node number
!     enddo
!
!     here we first setup the ia array
!
    ia(1) = 1
    k = 0
    do ipoin = 0, maxpatches
      if (ic(ipoin) .gt. 0) then ! skip empty colours
        k = k + 1
        if (k .gt. nclr) then
          write (6, *) 'setbndrynodeptr: found too many colours; check&
          & nclr !'
          error stop 1
        end if
        ia(k + 1) = ia(k) + ic(ipoin)
        iclr(k) = ipoin ! colour of the k-th boundary patch
      end if
    end do
    nnzr = ia(nclr + 1) - ia(1)
    write (6, *) 'setbndrynodeptr: ', nnzr, ' entries expected in ja'
!
!     allocate ja
!
    allocate (ja(nnzr), source=0_i4)
!
    write (6, *) 'setbndrynodeptr: finished initializing'
    write (6, *) 'setbndrynodeptr: now calling setbndrynodelist'
    write (6, *)
!
!     here we fill the ja array
!
    call setbndrynodelist(ia, ja, iclr, nclr, closed, bndfac, nbfac, nodptr, nbpoin)
  end subroutine setbndrynodeptr

end module mod_setbndrynodeptr
