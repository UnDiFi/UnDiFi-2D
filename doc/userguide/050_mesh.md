\hypertarget{mesh}{}

# Mesh Software \label{chap:mesh}
<!-- In this chapter, the description of the mesh generation process
for UnDiFi is detailed. -->

In the present version of the UnDiFi code, the
[Triangle](https://www.cs.cmu.edu/~quake/triangle.html) mesh generator
has been choosen. Since detailed information about Triangle can be found in the given link,
and since no changes were done to the original code, we point the
reader to the main website.

**Note:** future releases of the code may include several other mesh
generators such as [Tetgen](http://wias-berlin.de/software/tetgen/),
[Yams](https://www.ljll.math.upmc.fr/frey/software.html),
[Delaundo](http://www.ae.metu.edu.tr/tuncer/ae546/prj/delaundo/).

<!-- In the following, only those features which are used in the mesh
generation process with UnDiFi are described. -->

<!--

## Mesh generation/conversion with Triangle
Triangle is a C program which generates meshes, Delaunay triangulations
and Voronoi diagrams for 2D pointsets, by Jonathan Shewchuk. It
generates exact Delaunay triangulations, constrained Delaunay
triangulations, and quality conforming Delaunay triangulations. The
latter can be generated with no small angles, and are thus suitable for
finite element analysis. It produces its own documentation. Complete
instructions are printed by invoking the program with the `-h` switch:

        triangle -h

The instructions are long; you'll probably want to pipe the output to more or redirect it to a file. The programs gives a short
list of command line options if it is invoked without arguments You have to give TRIANGLE some input to work with. Usually, you
start with a Planar Straight Line Graph (PSLG) which defines a set of points, and some nonintersecting segments that you want to
be sure are included in the triangulation. This information can easily be stored in a POLY, which is one of the expected input
formats Normally, TRIANGLE will triangulate the entire region defined by the convex hull of the points. You may define a portion
of the convex hull to be a `hole`, which is not to be triangulated. This is done by adding some hole information to your POLY
file

Try out TRIANGLE on the sample file, A.poly:

        triangle -p A.poly

Triangle will read the Planar Straight Line Graph defined by A.poly, and write its constrained Delaunay triangulation to A.1.node,
A.1.ele and A.1.poly. A is a planar straight line graph of the capital letter A. We use it as input to get a constrained Delaunay
triangulation.

    A.poly, the POLY file;
    A.1.poly, the output polygon file;
    A.1.node, the output node file;
    A.1.ele; the output element file;

For contrast, try running

        triangle -pq A.poly

Now, click on the same `ele` button. A new triangulation will be loaded; this one having no angles smaller than 20 degrees.

To see a Voronoi diagram, try this:

        cp A.poly A.node
        triangle -v A

### .node files

First line: `<# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>`

Remaining lines: `<vertex #> <x> <y> [attributes] [boundary marker]`

Blank lines and comments prefixed by `#` may be placed anywhere. Vertices must be numbered consecutively, starting from one or
zero.

The attributes, which are typically floating-point values of physical quantities (such as mass or conductivity) associated with
the nodes of a finite element mesh, are copied unchanged to the output mesh. If -q, -a, -u, or -s is selected, each new Steiner
point added to the mesh will have quantities assigned to it by linear interpolation.

If the fourth entry of the first line is `1`, the last column of the remainder of the file is assumed to contain boundary markers.
Boundary markers are used to identify boundary vertices and vertices resting on PSLG segments. The .node files produced by
Triangle contain boundary markers in the last column unless they are suppressed by the -B switch.

### .ele files

First line: `<# of triangles> <nodes per triangle> <# of attributes>`
Remaining lines: `<triangle #> <node> <node> <node> ... [attributes]`

Blank lines and comments prefixed by `#` may be placed anywhere. Triangles must be numbered consecutively, starting from one
or zero. Nodes are indices into the corresponding .node file. The first three nodes are the corner vertices, and are listed in
counterclockwise order around each triangle. (The remaining nodes, if any, depend on the type of finite element used.)

As in .node files, the attributes are typically floating-point values of physical quantities (such as mass or conductivity)
associated with the elements (triangles) of a finite element mesh. Because there is no simple mapping from input to output
triangles, an attempt is made to interpolate attributes, which may result in a good deal of diffusion of attributes among
nearby triangles as the triangulation is refined. Attributes do not diffuse across segments, so attributes used to identify
segment-bounded regions remain intact.

In output .ele files, all triangles have three nodes each unless the -o2 switch is used, in which case subparametric quadratic
elements with six nodes are generated. The fourth, fifth, and sixth nodes lie on the midpoints of the edges opposite the first,
second, and third vertices.


### .poly files

    First line: <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>
    Following lines: <vertex #> <x> <y> [attributes] [boundary marker]
    One line: <# of segments> <# of boundary markers (0 or 1)>
    Following lines: <segment #> <endpoint> <endpoint> [boundary marker]
    One line: <# of holes>
    Following lines: <hole #> <x> <y>
    Optional line: <# of regional attributes and/or area constraints>
    Optional following lines: <region #> <x> <y> <attribute> <maximum area>

A .poly file represents a PSLG, as well as some additional information. PSLG stands for Planar Straight Line Graph, a term
familiar to computational geometers. By definition, a PSLG is just a list of vertices and segments. A .poly file can also contain
information about holes and concavities, as well as regional attributes and constraints on the areas of triangles.

The first section lists all the vertices, and is identical to the format of .node files. <# of vertices> may be set to zero to
indicate that the vertices are listed in a separate .node file; .poly files produced by Triangle always have this format. A vertex
set represented this way has the advantage that it may easily be triangulated with or without segments (depending on whether the
.poly or .node file is read).

The second section lists the segments. Segments are edges whose presence in the triangulation is enforced (although each segment
may be subdivided into smaller edges). Each segment is specified by listing the indices of its two endpoints. This means that you
must include its endpoints in the vertex list. Each segment, like each vertex, may have a boundary marker.

If -q, -a, -u, and -s are not selected, Triangle will produce a constrained Delaunay triangulation (CDT), in which each segment
appears as a single edge in the triangulation. If -q, -a, -u, or -s is selected, Triangle will produce a conforming constrained
Delaunay triangulation, in which segments may be subdivided into smaller edges. If -D is selected as well, Triangle will produce a
conforming Delaunay triangulation, so every triangle is Delaunay, and not just constrained Delaunay.

The third section lists holes (and concavities, if -c is selected) in the triangulation. Holes are specified by identifying a
point inside each hole. After the triangulation is formed, Triangle creates holes by eating triangles, spreading out from each
hole point until its progress is blocked by PSLG segments; you must be careful to enclose each hole in segments, or your whole
triangulation might be eaten away. If the two triangles abutting a segment are eaten, the segment itself is also eaten. Do not
place a hole directly on a segment; if you do, Triangle will choose one side of the segment arbitrarily.

The optional fourth section lists regional attributes (to be assigned to all triangles in a region) and regional constraints
on the maximum triangle area. Triangle will read this section only if the -A switch is used or the -a switch is used without a
number following it, and the -r switch is not used. Regional attributes and area constraints are propagated in the same manner
as holes; you specify a point for each attribute and/or constraint, and the attribute and/or constraint will affect the whole
region (bounded by segments) containing the point. If two values are written on a line after the x and y coordinates, the first
such value is assumed to be a regional attribute (but will only be applied if the -A switch is selected), and the second value
is assumed to be a regional area constraint (but will only be applied if the -a switch is selected). You may specify just one
value after the coordinates, which can serve as both an attribute and an area constraint, depending on the choice of switches.
If you are using the -A and -a switches simultaneously and wish to assign an attribute to some region without imposing an area
constraint, use a negative maximum area.

Blank lines and comments prefixed by `#` may be placed anywhere. Vertices, segments, holes, and regions must be numbered
consecutively, starting from one or zero. (The choice to begin the numbering from one or zero must be consistent across all
objects.)

When a triangulation is created from a .poly file, you must either enclose the entire region to be triangulated in PSLG segments,
or use the -c switch, which encloses the convex hull of the input vertex set. If you do not use the -c switch, Triangle will eat
all triangles that are not enclosed by segments; if you are not careful, your whole triangulation may be eaten away. If you do use
the -c switch, you can still produce concavities by the appropriate placement of holes just within the convex hull.

An ideal PSLG has no intersecting segments, nor any vertices that lie upon segments (except, of course, the endpoints of each
segment.) You aren't required to make your .poly files ideal, but you should be aware of what can go wrong. Segment intersections
are relatively safe - Triangle will calculate the intersection points for you and add them to the triangulation - as long as
your machine's floating-point precision doesn't become a problem. You are tempting the fates if you have three segments that
cross at the same location, and expect Triangle to figure out where the intersection point is. Thanks to floating-point roundoff
error, Triangle will probably decide that the three segments intersect at three different points, and you will find a minuscule
triangle in your output - unless Triangle tries to refine the tiny triangle, uses up the last bit of machine precision, and fails
to terminate at all. You're better off putting the intersection point in the input files, and manually breaking up each segment
into two. Similarly, if you place a vertex at the middle of a segment, and hope that Triangle will break up the segment at that
vertex, you might get lucky. On the other hand, Triangle might decide that the vertex doesn't lie precisely on the line, and
you'll have a needle-sharp triangle in your output - or a lot of tiny triangles if you're generating a quality mesh.

When Triangle reads a .poly file, it also writes a .poly file, which includes all edges that are subsegments of input segments.
If the -c switch is used, the output .poly file will also include all of the edges on the convex hull. Hence, the output .poly
file is useful for finding edges associated with input segments and setting boundary conditions in finite element simulations.
More importantly, you will need it if you plan to refine the output mesh, and don't want segments to be missing in later
triangulations.

### .edge files

First line: `<# of edges> <# of boundary markers (0 or 1)>`
Following lines: `<edge #> <endpoint> <endpoint> [boundary marker]`

Blank lines and comments prefixed by `#` may be placed anywhere. Edges are numbered consecutively, starting from one or zero.
Endpoints are indices into the corresponding .node file.

Triangle can produce .edge files (use the -e switch), but cannot read them. The optional column of boundary markers is suppressed
by the -B switch.

In Voronoi diagrams, one also finds a special kind of edge that is an infinite ray with only one endpoint. For these edges, a
different format is used:

    <edge #> <endpoint> -1 <direction x> <direction y>

The `direction` is a floating-point vector that indicates the direction of the infinite ray.

### .neigh files

First line: `<# of triangles> <# of neighbors per triangle (always 3)>`
Following lines: `<triangle #> <neighbor> <neighbor> <neighbor>`

Blank lines and comments prefixed by `#` may be placed anywhere. Triangles are numbered consecutively, starting from one or zero.
Neighbors are indices into the corresponding .ele file. An index of -1 indicates no neighbor (because the triangle is on an
exterior boundary). The first neighbor of triangle i is opposite the first corner of triangle i, and so on.

Triangle can produce .neigh files (use the -n switch), but cannot read them.

-->

<!--

### Command line switches

To run Triangle, the command line syntax is

`triangle [-prq__a__uAcDjevngBPNEIOXzo_YS__iFlsCQVh] input_file`

As in .node files, the attributes are typically floating-point values of physical quantities (such as mass or conductivity)
associated with the elements (triangles) of a finite element mesh. Because there is no simple mapping from input to output
triangles, an attempt is made to interpolate attributes, which may result in a good deal of diffusion of attributes among
nearby triangles as the triangulation is refined. Attributes do not diffuse across segments, so attributes used to identify
segment-bounded regions remain intact.

```
-p Triangulates a Planar Straight Line Graph (.poly file).
-r Refines a previously generated mesh.
-q Quality mesh generation with no angles smaller than 20 degrees. An alternate minimum angle may be specified after the `q'.
-a Imposes a maximum triangle area constraint. A fixed area constraint (that applies to every triangle) may be specified after the `a', or varying area constraints may be read from a .poly file or .area file.
-u Imposes a user-defined constraint on triangle size.
-A Assigns a regional attribute to each triangle that identifies what segment-bounded region it belongs to.
-c Encloses the convex hull with segments.
-D Conforming Delaunay: use this switch if you want all triangles in the mesh to be Delaunay, and not just constrained Delaunay; or if you want to ensure that all Voronoi vertices lie within the triangulation.
-j Jettisons vertices that are not part of the final triangulation from the output .node file (including duplicate input vertices and vertices ``eaten'' by holes).
-e Outputs (to an .edge file) a list of edges of the triangulation.
-v Outputs the Voronoi diagram associated with the triangulation. Does not attempt to detect degeneracies, so some Voronoi vertices may be duplicated.
-n Outputs (to a .neigh file) a list of triangles neighboring each triangle.
-g Outputs the mesh to an Object File Format (.off) file, suitable for viewing with the Geometry Center's Geomview package.
-B Suppresses boundary markers in the output .node, .poly, and .edge output files.
-P Suppresses the output .poly file. Saves disk space, but you lose the ability to maintain constraining segments on later refinements of the mesh.
-N Suppresses the output .node file.
-E Suppresses the output .ele file.
-I Suppresses mesh iteration numbers.
-O Suppresses holes: ignores the holes in the .poly file.
-X Suppresses exact arithmetic.
-z Numbers all items starting from zero (rather than one). Note that this switch is normally overrided by the value used to number the first vertex of the input .node or .poly file. However, this switch is useful when calling Triangle from another program.
-o2 Generates second-order subparametric elements with six nodes each.
-Y Prohibits the insertion of Steiner points on the mesh boundary. If specified twice (-YY), it prohibits the insertion of Steiner points on any segment, including internal segments.
-S Specifies the maximum number of added Steiner points.
-i Uses the incremental algorithm for Delaunay triangulation, rather than the divide-and-conquer algorithm.
-F Uses Steven Fortune's sweepline algorithm for Delaunay triangulation, rather than the divide-and-conquer algorithm.
-l Uses only vertical cuts in the divide-and-conquer algorithm. By default, Triangle uses alternating vertical and horizontal cuts, which usually improve the speed except with vertex sets that are small or short and wide. This switch is primarily of theoretical interest.
-s Specifies that segments should be forced into the triangulation by recursively splitting them at their midpoints, rather than by generating a constrained Delaunay triangulation. Segment splitting is true to Ruppert's original algorithm, but can create needlessly small triangles. This switch is primarily of theoretical interest.
-C Check the consistency of the final mesh. Uses exact arithmetic for checking, even if the -X switch is used. Useful if you suspect Triangle is buggy.
-Q Quiet: Suppresses all explanation of what Triangle is doing, unless an error occurs.
-V Verbose: Gives detailed information about what Triangle is doing. Add more `V's for increasing amount of detail. `-V' gives information on algorithmic progress and detailed statistics.
-h Help: Displays complete instructions.
```

### POLY
POLY is a data directory which contains examples of POLY files, a
format used by Jonathan Shewchuk to define PSLG's, planar straight line
graphs, for use with his program TRIANGLE. A Planar Straight Line Graph
(PSLG) is a set of vertices and segments. Segments are simply edges,
whose endpoints are vertices in the PSLG. Segments may intersect each
other only at their endpoints.

POLY File Characteristics:

  * ASCII
  * 2D
  * vertices are specified by coordinates.
  * line segments are specified by listing the indices of pairs of vertices.
  * a hole may be specified by listing the coordinates of a point inside the hole.
  * No compression
  * 1 image

Comments are prefixed by the character `#`. Everything from the comment character to the end of the line is ignored.
Vertices, segments, holes, and regions must be numbered and listed consecutively, starting from either 1 or 0.

The first line lists

  * The number of vertices (this is sometimes set to 0, to indicate that the vertices should be read from a NODE file);
  * The spatial dimension, which must be 2;
  * The number of vertex attributes;
  * The number of vertex boundary markers, which must be 0 or 1.

The vertex records must follow, with the format:

  * vertex index (these must be consecutive, starting either from 0 or 1);
  * X and Y coordinates;
  * The vertex attributes (if any);
  * The vertex boundary marker (if any).

The next line lists

  * The number of segments;
  * The number of segment boundary markers (0 or 1).

Segments should not cross each other; vertices should only lie on the ends of segments, and are never contained inside a segment.

The segments records must follow, with the format:

  * segment index;
  * start vertex, end vertex;
  * Boundary marker (if any).

The third section lists holes (and concavities, if -c is selected) in
the triangulation. Holes are specified by identifying a point inside
each hole. After the triangulation is formed, Triangle creates holes by
eating triangles, spreading out from each hole point until its progress
is blocked by PSLG segments; you must be careful to enclose each hole
in segments, or your whole triangulation might be eaten away. If the
two triangles abutting a segment are eaten, the segment itself is also
eaten. Do not place a hole directly on a segment; if you do, Triangle
chooses one side of the segment arbitrarily.

The next line lists

  * The number of holes.

The hole records must follow, with the format:

  * hole index;
  * X coordinate, Y coordinate of some point within the hole.

The optional fourth section lists regional attributes (to be assigned
to all triangles in a region) and regional constraints on the maximum
triangle area. Triangle reads this section only if the -A switch is
used or the -a switch is used without a number following it, and the
-r switch is not used. Regional attributes and area constraints are
propagated in the same manner as holes; you specify a point for each
attribute and/or constraint, and the attribute and/or constraint
affects the whole region (bounded by segments) containing the point.
If two values are written on a line after the x and y coordinate, the
first such value is assumed to be a regional attribute (but is only
applied if the -A switch is selected), and the second value is assumed
to be a regional area constraint (but is only applied if the -a switch
is selected). You may specify just one value after the coordinates,
which can serve as both an attribute and an area constraint, depending
on the choice of switches. If you are using the -A and -a switches
simultaneously and wish to assign an attribute to some region without
imposing an area constraint, use a negative maximum area.

The next line is optional. If given, it lists

  * The number of region attributes.

The optional regional attributes records must follow, with the format:

  * region index;
  * X coordinate, Y coordinate of a point in the region;
  * Attributes (if any);
  * Maximum area of triangles in the region;

A Sample POLY file:

Here is a sample file box.poly describing a square with a square hole.

~~~~~~
# A box with eight vertices in 2D, no attributes, one boundary marker.
8 2 0 1
# Outer box has these vertices:
 1   0 0   0
 2   0 3   0
 3   3 0   0
 4   3 3   33     # A special marker for this vertex.
# Inner square has these vertices:
 5   1 1   0
 6   1 2   0
 7   2 1   0
 8   2 2   0
# Five segments with boundary markers.
5 1
 1   1 2   5      # Left side of outer box.
# Square hole has these segments:
 2   5 7   0
 3   7 8   0
 4   8 6   10
 5   6 5   0
# One hole in the middle of the inner square.
1
1   1.5 1.5
~~~~~~

<!-- ## Mesh generation/conversion with Tetgen   -->
<!-- ## Mesh generation/conversion with Yams     -->
<!-- ## Mesh generation/conversion with Delaundo -->
