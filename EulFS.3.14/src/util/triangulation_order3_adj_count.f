      subroutine triangulation_order3_adj_count ( node_num, 
     &  triangle_num, triangle_node, triangle_neighbor, adj_num, 
     &  adj_col )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_ADJ_COUNT counts adjacencies in a triangulation.
c
c  Discussion:
c
c    This routine is called to count the adjacencies, so that the
c    appropriate amount of memory can be set aside for storage when
c    the adjacency structure is created.
c
c    The triangulation is assumed to involve 3-node triangles.
c
c    Two nodes are "adjacent" if they are both nodes in some triangle.
c    Also, a node is considered to be adjacent to itself.
c
c  Diagram:
c
c       3
c    s  |\
c    i  | \
c    d  |  \
c    e  |   \  side 2
c       |    \
c    3  |     \
c       |      \
c       1-------2
c
c         side 1
c
c    The local node numbering
c
c
c   21-22-23-24-25
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c   16-17-18-19-20
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c   11-12-13-14-15
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c    6--7--8--9-10
c    |\ |\ |\ |\ |
c    | \| \| \| \|
c    1--2--3--4--5
c
c    A sample grid.
c
c
c    Below, we have a chart that summarizes the adjacency relationships
c    in the sample grid.  On the left, we list the node, and its neighbors,
c    with an asterisk to indicate the adjacency of the node to itself
c    (in some cases, you want to count this self adjacency and in some
c    you don't).  On the right, we list the number of adjacencies to
c    lower-indexed nodes, to the node itself, to higher-indexed nodes,
c    the total number of adjacencies for this node, and the location
c    of the first and last entries required to list this set of adjacencies
c    in a single list of all the adjacencies.
c
c    N   Adjacencies                Below  Self   Above   Total First  Last
c
c   --  -- -- -- -- -- -- --           --    --      --      --   ---     0   
c    1:  *  2  6                        0     1       2       3     1     3
c    2:  1  *  3  6  7                  1     1       3       5     4     8
c    3:  2  *  4  7  8                  1     1       3       5     9    13
c    4:  3  *  5  8  9                  1     1       3       5    14    18
c    5:  4  *  9 10                     1     1       2       4    19    22
c    6:  1  2  *  7 11                  2     1       2       5    23    27
c    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
c    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
c    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
c   10:  5  9  * 14 15                  2     1       2       5    49    53
c   11:  6  7  * 12 16                  2     1       2       5    54    58
c   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
c   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
c   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
c   15: 10 14  * 19 20                  2     1       2       5    80    84
c   16: 11 12  * 17 21                  2     1       2       5    85    89
c   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
c   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
c   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
c   20: 15 19  * 24 25                  2     1       2       5   111   115
c   21: 16 17  * 22                     2     1       1       4   116   119
c   22: 17 18 21  * 23                  3     1       1       5   120   124
c   23: 18 19 22  * 24                  3     1       1       5   125   129
c   24: 19 20 23  * 25                  3     1       1       5   130   134
c   25: 20 24  *                        2     1       0       3   135   137
c   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
c      
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(3,TRIANGLE_NUM), lists the 
c    nodes that make up each triangle, in counterclockwise order. 
c
c    Input, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), for each 
c    side of a triangle, lists the neighboring triangle, or -1 if there is
c    no neighbor.
c
c    Output, integer ADJ_NUM, the number of adjacencies.
c
c    Output, integer ADJ_COL(NODE_NUM+1).  Information about 
c    column J is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
c
      implicit none

      integer node_num
      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer adj_num
      integer adj_col(node_num+1)
      integer i
      integer n1
      integer n2
      integer n3
      integer triangle
      integer triangle2
      integer triangle_neighbor(3,triangle_num)
      integer triangle_node(triangle_order,triangle_num)

      adj_num = 0
c
c  Set every node to be adjacent to itself.
c
      do i = 1, node_num
        adj_col(i) = 1
      end do
caldo
      do triangle = 1, triangle_num
         do i = 1,triangle_order
            triangle2 = triangle_neighbor(i,triangle)
            if( triangle2 .gt. triangle_num )then
                triangle_neighbor(i,triangle) = -triangle2
            endif
         enddo
      enddo
caldo
c
c  Examine each triangle.
c
      do triangle = 1, triangle_num

        n1 = triangle_node(1,triangle)
        n2 = triangle_node(2,triangle)
        n3 = triangle_node(3,triangle)
c
c  Add edge (1,2) if this is the first occurrence,
c  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
c  or if this triangle is the first of the pair in which the edge
c  occurs (TRIANGLE < TRIANGLE2).
c
        triangle2 = triangle_neighbor(1,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then
          adj_col(n1) = adj_col(n1) + 1
          adj_col(n2) = adj_col(n2) + 1
        end if
c
c  Add edge (2,3).
c
        triangle2 = triangle_neighbor(2,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then
          adj_col(n2) = adj_col(n2) + 1
          adj_col(n3) = adj_col(n3) + 1
        end if
c
c  Add edge (3,1).
c
        triangle2 = triangle_neighbor(3,triangle)

        if ( triangle2 .lt. 0 .or. triangle .lt. triangle2 ) then
          adj_col(n1) = adj_col(n1) + 1
          adj_col(n3) = adj_col(n3) + 1
        end if
          
      end do
c
c  We used ADJ_COL to count the number of entries in each column.
c  Convert it to pointers into the ADJ array.
c
      do i = node_num, 1, -1
        adj_col(i+1) = adj_col(i)
      end do

      adj_col(1) = 1
      do i = 2, node_num+1
        adj_col(i) = adj_col(i-1) + adj_col(i)
      end do

      adj_num = adj_col(node_num+1) - 1

      return
      end
