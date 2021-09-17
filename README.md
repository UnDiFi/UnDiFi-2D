UnDiFi-2D
=========

UnDiFi-2D: an Unstructured Discontinuity Fitting code for 2D grids.

It is written in standard (compliant) Fortran 77/95 with highly
modularity as design target.

The aim of UnDiFi-2D is to explicitely manage discontinuities in the flow
field. In our unstructured shock-fitting approach the shock front is
described using a double-sided, polygonal curve. Two sets of flow
states, corresponding to the upstream and downstream sides of the
discontinuity, are assigned to the grid-points located on either side
of the shock front. This is allowed to move, while obeying to the
Rankine-Hugoniot jump relations, throughout a background triangular
mesh that covers the entire computational domain. At each time step,
a local, constrained Delaunay triangulation is applied in the
neighbourhood of the shock front to ensure that the edges that make up
the shock front are part of the overall triangular grid. The fitted
shock acts as an interior boundary for the shock-capturing solver that
is used to solve the discretised governing equations in the smooth
regions of the flow-field.

Copyrights
----------

UnDiFi-2D is an open source project, it is distributed under the
[GPL v3](http://www.gnu.org/licenses/gpl-3.0.html). Anyone is interest
to use, to develop or to contribute to UnDiFi-2D is welcome.

Documentation
-------------

Detailed documentation can be found in the `doc` folder.

Citing
------

If you publish work which mentions UnDiFi-2D, or UnDiFi-2D has been useful
in your research, please kindly cite the following paper:

~~~
@article{campoli-2019,
  author  = {L. Campoli, A. Assonitis, A. Bonfiglioli, Ciallella, R. Paciorri, M. Ricchiuto},
  title   = {{UnDiFi, an Unstructured Discontinuity Fitting algorithm code}},
  journal = {Computer Physics Communications },
  volume  = {},
  number  = {},
  pages   = {},
  year    = {},
  issn    = {},
  doi     = {},
  url     = {},
}
~~~

Bibliography
------------

[1] Aldo Bonfiglioli. Fluctuation splitting schemes for the compressible
and incompressible euler and navier-stokes equations. International
Journal of Computational Fluid Dynamics, 14(1):21–39, 2000.

[2] Andrea Lani, Tiago Quintino, Dries Kimpe, Herman Deconinck, Stefan
Vandewalle, and Stefaan Poedts. The coolfluid framework: design
solutions for high performance object oriented scientific computing
software. In International Conference on Computational Science, pages
279–286. Springer, 2005.

[3] Mario Ricchiuto and Remi Abgrall. Explicit runge–kutta residual
distribution schemes for time dependent problems: second order case.
Journal of Computational Physics, 229(16):5653–5691, 2010.

[4] H Deconinck, H Paillere, R Struijs, and Philip L Roe.
Multidimensional upwind schemes based on fluctuation-splitting
for systems of conservation laws. Computational Mechanics,
11(5-6):323–340, 1993.
