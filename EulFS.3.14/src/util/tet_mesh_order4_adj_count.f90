subroutine tet_mesh_order4_adj_count ( node_num, tetra_num, tetra_node, &
  adj_num, adj_row )

!*****************************************************************************80
!
!! TET_MESH_ORDER4_ADJ_COUNT counts the number of nodal adjacencies.
!
!  Discussion:
!
!    Assuming that the tet mesh is to be used in a finite element
!    computation, we declare that two distinct nodes are "adjacent" if and
!    only if they are both included in some tetrahedron.
!
!    It is the purpose of this routine to determine the number of
!    such adjacency relationships.
!
!    The initial count gets only the (I,J) relationships, for which
!    node I is strictly less than node J.  This value is doubled
!    to account for symmetry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) TETRA_NUM, the number of tetrahedrons.
!
!    Input, integer ( kind = 4 ) TETRA_NODE(4,TETRA_NUM), the indices of the nodes.
!
!    Output, integer ( kind = 4 ) ADJ_NUM, the total number of adjacency relationships,
!
!    Output, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1), the ADJ pointer array.
!
  implicit none

  integer ( kind = 4 ) tetra_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) adj_row(node_num+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) pair(2,6*tetra_num)
  integer ( kind = 4 ) pair_num
  integer ( kind = 4 ) pair_unique_num
  integer ( kind = 4 ) tetra_node(4,tetra_num)
!
!  Each order 4 tetrahedron defines 6 adjacency pairs.
!
  pair(1,            1:  tetra_num) = tetra_node(1,1:tetra_num)
  pair(2,            1:  tetra_num) = tetra_node(2,1:tetra_num)

  pair(1,  tetra_num+1:2*tetra_num) = tetra_node(1,1:tetra_num)
  pair(2,  tetra_num+1:2*tetra_num) = tetra_node(3,1:tetra_num)

  pair(1,2*tetra_num+1:3*tetra_num) = tetra_node(1,1:tetra_num)
  pair(2,2*tetra_num+1:3*tetra_num) = tetra_node(4,1:tetra_num)

  pair(1,3*tetra_num+1:4*tetra_num) = tetra_node(2,1:tetra_num)
  pair(2,3*tetra_num+1:4*tetra_num) = tetra_node(3,1:tetra_num)

  pair(1,4*tetra_num+1:5*tetra_num) = tetra_node(2,1:tetra_num)
  pair(2,4*tetra_num+1:5*tetra_num) = tetra_node(4,1:tetra_num)

  pair(1,5*tetra_num+1:6*tetra_num) = tetra_node(3,1:tetra_num)
  pair(2,5*tetra_num+1:6*tetra_num) = tetra_node(4,1:tetra_num)

  pair_num = 6 * tetra_num
!
!  Force the nodes of each pair to be listed in ascending order.
!
  call i4col_sort2_a ( 2, pair_num, pair )
!
!  Rearrange the columns in ascending order.
!
  call i4col_sort_a ( 2, pair_num, pair )
!
!  Get the number of unique columns.
!
  call i4col_sorted_unique_count ( 2, pair_num, pair, pair_unique_num )
!
!  The number of adjacencies is TWICE this value, plus the number of nodes.
!
  adj_num = 2 * pair_unique_num
!
!  Now set up the ADJ_ROW counts.
!
  adj_row(1:node_num) = 0

  do k = 1, pair_num

    if ( 1 < k ) then
      if ( pair(1,k-1) == pair(1,k) .and. &
           pair(2,k-1) == pair(2,k) ) then
        cycle
      end if
    end if

    i = pair(1,k)
    j = pair(2,k)

    adj_row(i) = adj_row(i) + 1
    adj_row(j) = adj_row(j) + 1

  end do
!
!  We used ADJ_ROW to count the number of entries in each row.
!  Convert it to pointers into the ADJ array.
!
  adj_row(2:node_num+1) = adj_row(1:node_num)

  adj_row(1) = 1
  do i = 2, node_num+1
    adj_row(i) = adj_row(i-1) + adj_row(i)
  end do

  return
end
subroutine i4col_sort_a ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT_A ascending sorts an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of integer values, regarded
!    as an array of N columns of length M.
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in ascending
!    lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4col_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4col_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4col_sort2_a ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT2_A ascending sorts the elements of each column of an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of integer values, regarded
!    as an array of N columns of length M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A, and the length
!    of a vector of data.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M vectors.
!    On output, the elements of each column of A have been sorted in ascending
!    order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 1 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if
!
!  Initialize.
!
  do col = 1, n

    i = 0
    indx = 0
    isgn = 0
    j = 0
!
!  Call the external heap sorter.
!
    do

      call sort_heap_external ( m, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
      if ( 0 < indx ) then

        call i4_swap ( a(i,col), a(j,col) )
!
!  Compare the I and J objects.
!
      else if ( indx < 0 ) then

        if ( a(j,col) < a(i,col) ) then
          isgn = +1
        else
          isgn = -1
        end if

      else if ( indx == 0 ) then

        exit

      end if

    end do

  end do

  return
end
subroutine i4col_sorted_unique_count ( m, n, a, unique_num )

!*****************************************************************************80
!
!! I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of integer values, regarded
!    as an array of N columns of length M.
!
!    The columns of the array may be ascending or descending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), a sorted array, containing
!    N columns of data.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1
  j1 = 1

  do j2 = 2, n

    if ( any ( a(1:m,j1) /= a(1:m,j2) ) ) then
      unique_num = unique_num + 1
      j1 = j2
    end if

  end do

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
subroutine i4col_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4COL_COMPARE compares columns I and J of an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of integer values, regarded
!    as an array of N columns of length M.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of N columns of vectors of length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column J < column I.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Check.
!
  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a,i8,a)' ) '  Column index I = ', i, ' is less than 1.'
    stop
  end if

  if ( n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a,i8,a)' ) '  N = ', n, ' is less than column index I = ', i
    stop
  end if

  if ( j < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a,i8,a)' ) '  Column index J = ', j, ' is less than 1.'
    stop
  end if

  if ( n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a,i8,a)' ) '  N = ', n, ' is less than column index J = ', j
    stop
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      isgn = -1
      return
    else if ( a(k,j) < a(k,i) ) then
      isgn = +1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements I and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

  return
end
subroutine i4col_swap ( m, n, a, j1, j2 )

!*****************************************************************************80
!
!! I4COL_SWAP swaps columns J1 and J2 of an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of integer values, regarded
!    as an array of N columns of length M.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, J1 = 2, J2 = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      A = (
!        1  4  3  2
!        5  8  7  6
!        9 12 11 10 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in the array.
!
!    Input/output, integer ( kind = 4 ) A(M,N), an array of N columns of length M.
!
!    Input, integer ( kind = 4 ) J1, J2, the columns to be swapped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col(m)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2

  if ( j1 < 1 .or. n < j1 .or. j2 < 1 .or. n < j2 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  J1 or J2 is out of bounds.'
    write ( *, '(a,i8)' ) '  J1 =    ', j1
    write ( *, '(a,i8)' ) '  J2 =    ', j2
    write ( *, '(a,i8)' ) '  N =     ', n
    stop

  end if

  if ( j1 == j2 ) then
    return
  end if

  col(1:m)  = a(1:m,j1)
  a(1:m,j1) = a(1:m,j2)
  a(1:m,j2) = col(1:m)

  return
end
