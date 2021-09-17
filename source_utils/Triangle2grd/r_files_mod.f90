module r_files_mod

! =====================================
! Thorough reading of different files
! and writing in arrays
! =====================================

  implicit none

  contains

    subroutine r_files(xy_node,att_node,ele,neigh,b_edge,bmk_f,nb_bed,fname,lenfname,type,bmk_edge)

      real*8,dimension(:,:),intent(inout):: xy_node, att_node
      integer,dimension(:,:),intent(inout):: b_edge, ele, neigh, bmk_f
      integer,dimension(:),intent(inout)::bmk_edge
      integer,dimension(:),intent(in):: lenfname, type

      integer,intent(in)::nb_bed

      character(len=*),dimension(:),intent(in)::fname

!aldo
      integer,dimension(:):: howmany(10)
!aldo
      integer::nb_node,nb_ele,nb_edge,nb_att, s1, s2, mark, nb_poly
      integer::i, j, temp, q, b ,ibfac, ielem

      do i = 1,10
         howmany(i) = 0
      enddo

      ! ======================
      ! Reading of .node file
      ! ======================
      open(unit=10, file=fname(1)(1:lenfname(1)), status='old')
      read(10,*) nb_node, temp, nb_att
      do i=1,nb_node
         read(10,*) temp, xy_node(1,i), xy_node(2,i), (att_node(j,i),j=1,nb_att)
      end do
      close(10)

      ! ==========================
      ! Reading of .ele node file
      ! ==========================
      open(11, file=fname(2)(1:lenfname(2)), status='old')
      read(11,*) nb_ele
      do i=1, nb_ele
         read(11,*) temp, ele(1,i), ele(2,i), ele(3,i)
      end do
      close(11)

      ! ============================
      ! Reading of .neigh node file
      ! ============================
      open(12, file=fname(3)(1:lenfname(3)), status='old')
      read(12,*)
      do i=1,nb_ele
         read(12,*) temp, neigh(1,i), neigh(2,i), neigh(3,i)
      end do
      close(12)

! *******  Nombre de faces d'element sans voisins (lorsqu'il y a -1) verification *********
!!$      IBFAC = 0
!!$      DO  IELEM = 1, nb_ele
!!$         DO  I = 1,3
!!$            IF( neigh(I,IELEM) .EQ. -1 )THEN
!!$               IBFAC = IBFAC+1
!!$            ENDIF
!!$         end DO
!!$      end DO
!!$!print*,'ibfac :', ibfac, 'nb_bedge',nb_bed
!!$if (ibfac/=nb_bed) then
!!$   print*,'problem ibfac/=nb_bed'
!!$   stop
!!$end if
! **********************

      ! ===========================
      ! Reading of .edge node file
      ! ===========================

      open(13, file=fname(4)(1:lenfname(4)), status='old')
      read(13,*) nb_edge
      j = 0
      do i=1, nb_edge
         read(13,*) temp, s1, s2, mark
         if ((mark/=0)) then
            j = j + 1                               
            b_edge(1,j) = s1
            b_edge(2,j) = s2
            bmk_edge(j) = mark
         end if
      end do
      close(13)

      ! ======================
      ! Boundary faces marker
      ! ======================
      
      do i = 1, nb_bed
         if (bmk_edge(i) == 1 )  bmk_f(:,i) = type(1)
         if (bmk_edge(i) == 2 )  bmk_f(:,i) = type(2)
         if (bmk_edge(i) == 3 )  bmk_f(:,i) = type(3)
         if (bmk_edge(i) == 4 )  bmk_f(:,i) = type(4)
         if (bmk_edge(i) == 5 )  bmk_f(:,i) = type(5)
         if (bmk_edge(i) == 6 )  bmk_f(:,i) = type(6)
         if (bmk_edge(i) == 7 )  bmk_f(:,i) = type(7)
         if (bmk_edge(i) == 8 )  bmk_f(:,i) = type(8)
         if (bmk_edge(i) == 9 )  bmk_f(:,i) = type(9)
         if (bmk_edge(i) == 10 ) bmk_f(:,i) = type(10)
         j = bmk_edge(i)
         howmany(j) = howmany(j) + 1
      end do
      do i = 1, 10
         j = howmany(i)
         if(j.NE.0)write(6,*)' subroutine r_files: faces coloured ',i,' are ',j
      end do

      do i = 1, nb_bed
         do j = i, nb_bed

            if (b_edge(1,i)==b_edge(1,j)) then
               if (bmk_f(1,i)/=bmk_f(1,j)) then
                  if ((bmk_f(1,i)==3) .or. bmk_f(1,j)==3) then
                     bmk_f(1,i) = 3
                     bmk_f(1,j) = 3
                  end if
                  if ((bmk_f(1,i)==1) .or. bmk_f(1,j)==1) then
                     bmk_f(1,i) = 1
                     bmk_f(1,j) = 1
                  end if
               end if
            end if

            if (b_edge(1,i)==b_edge(2,j)) then
               if (bmk_f(1,i)/=bmk_f(2,j)) then
                  if ((bmk_f(1,i)==3) .or. bmk_f(2,j)==3) then
                     bmk_f(1,i) = 3
                     bmk_f(2,j) = 3
                  end if
                  if ((bmk_f(1,i)==1) .or. bmk_f(2,j)==1) then
                     bmk_f(1,i) = 1
                     bmk_f(2,j) = 1
                  end if
               end if
            end if

            if (b_edge(2,i)==b_edge(1,j)) then
               if (bmk_f(2,i)/=bmk_f(1,j)) then
                  if ((bmk_f(2,i)==3) .or. bmk_f(1,j)==3) then
                     bmk_f(2,i) = 3
                     bmk_f(1,j) = 3
                  end if
                  if ((bmk_f(2,i)==1) .or. bmk_f(1,j)==1) then
                     bmk_f(2,i) = 1
                     bmk_f(1,j) = 1
                  end if
               end if
            end if

            if (b_edge(2,i)==b_edge(2,j)) then
               if (bmk_f(2,i)/=bmk_f(2,j)) then
                  if ((bmk_f(2,i)==3) .or. bmk_f(2,j)==3) then
                     bmk_f(2,i) = 3
                     bmk_f(2,j) = 3
                  end if
                  if ((bmk_f(2,i)==1) .or. bmk_f(2,j)==1) then
                     bmk_f(2,i) = 1
                     bmk_f(2,j) = 1
                  end if
               end if
            end if
         end do
      end do
      
    end subroutine r_files
    
  end module r_files_mod
  
