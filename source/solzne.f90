! Subroutine for reading(mode='r') and writing(mode='w') nodal values

subroutine solzne(filename, varray, nofvar, npoin, mode)

  use mod_kinds, only: wp, i4
  implicit none(type, external)

!     .. scalar arguments ..
  integer(i4) nofvar, npoin
  character filename*(*), mode*(*)

!     .. array arguments ..
  real(wp) varray(nofvar, npoin)

!     .. local scalars ..
  integer(i4) ifail, ixdrs, npold, nvold

!     .. external functions ..
  integer(i4) initxdr
  integer(i4) ixdrint, ixdrimat, ixdrclose, ixdrdmat
  external initxdr

!     .. external subroutines ..
  external seterr, ixdrint, ixdrimat, ixdrclose, ixdrdmat

!     .. data statements ..
  data ifail/0/

!     .. reading or backing up ..
!      if (mode.eq.'r') then
!        write (6,fmt=110) filename
!      else
!        write (6,fmt=112) filename
!      endif

  ixdrs = initxdr(filename, mode, .false.)

  npold = npoin
  nvold = nofvar
  ifail = ixdrint(ixdrs, npoin)
  ifail = ixdrint(ixdrs, nofvar)
  if (mode .eq. 'r') then
    if (npoin .ne. npold) call seterr&
    &(30hinconsistent npoin in datafile,&
    &30, 1, 2)
    if (nofvar .ne. nvold) call seterr&
    &(31hinconsistent nofvar in datafile,&
    &31, 1, 2)
  end if

  ifail = ixdrdmat(ixdrs, nofvar*npoin, varray)
!     call x04caf('general',' ',nofvar,npoin,varray,nofvar,
!    +            'nodal values',ifail)
  ifail = ixdrclose(ixdrs)

  return

110 format(/, 5x, 'reading solution from ', a40,/)
112 format(/, 5x, 'writing solution to ', a40,/)
115 format(' done',/)

end subroutine solzne
