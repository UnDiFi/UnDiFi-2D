\hypertarget{algorithm}{}

# Discontinuity Fitting Algorithm \label{chap:algorithm}

In this chapter, each step of the complete discontinuity fitting algorithm of **UnDiFi** is detailed.

## Steps of the Discontinuity Fitting Algorithm

The shock-fitting module identifies the grid-points and constrained edges (triangles in 3D) that make up the shock-fitting mesh
to be used over the current time interval and calls the meshing software that builds the shock-fitting mesh using a Constrained
Delauney Triangulation (CDT) (The meshing software).

Two different Fortran versions have been developed by the authors. The original F77 code is capable of dealing with steady
and unsteady two-dimensional flows, featuring both shocks and contact discontinuities as well as their mutual interactions or
their interaction with solid surfaces. A different F90 code is capable of dealing with steady three-dimensional flows featuring
multiple, but non-interacting shocks; it has been used to compute the results presented in [4]. A C++ implementation of the
shock-fitting module, capable of handling both two- and three-dimensional flows, has been developed at the Von Karman Institute;
it is described in [10] and it has also been made publicly accessible at: https://github.com/andrealani/ShockFitting

Steps of the Algorithm

The approach is inherently time-dependent: both the solution and the grid change with time, due to the displacement of the fitted
discontinuities. When a steady solution exists, the shock speed will asymptotically vanish and the tessellation of the flow domain
will not any longer change.

At time t the set of dependent variables and grid velocity are available within all grid-points of a tessellation (made of
triangles in 2D and tetrahedra in 3D) that covers the entire computational domain; this is what we call the "background" mesh.
In addition to the background mesh, the fitted discontinuities (either shocks or contact discontinuities) are discretized using
a collection of points which are mutually joined to form a connected series of line segments. Each fitted discontinuity is a
double-sided internal boundary of zero thickness. Because the width of the discontinuity is negligible, its two sides are
discretized using the same polygonal curve or triangulated surface; each pair of nodes that face each other on the two sides of
the discontinuity share the same geometrical location, but store different values of the dependent variables, one corresponding
to the upstream state and the other to the downstream one. Moreover, a velocity vector normal to the discontinuity is assigned to
each pair of grid-points on the fitted discontinuity, representing its displacement velocity. The spatial location of the fitted
discontinuities is independent of the location of the grid-points that make up the background grid.

The process that leads from the available mesh and solution at time t to an updated mesh and solution at time t+dt can be split
into seven steps that will be described in the following paragraphs and shown in Fig. undifi_algorithm_steps.

![Fig.1](./images/f1.png?raw=true)

Figure 1: UnDiFi Unstructured shock-fitting. (a) Shock front moving over the background triangular mesh at time t.
(b) Dashed lines mark the cells to be removed; dashed circles denote the phantom nodes.
(c) The background mesh is split into disjoint sub-domains by a hole which encloses the shock.
(d) The triangulation around the shock has been rebuilt. (e) Calculation of the shock-tangent and shock-normal unit vectors.
(f) The shock displacement induces mesh deformation.

## Cell Remoaval Around the Shock Front

In this first step, the fitted discontinuities are laid on top of the background mesh, as shown in Fig. undifi_algorithm_steps (a). All those cells that are crossed by the fitted discontinuities and those mesh points that are located too close to it are temporarily removed from the background mesh, as shown in Fig. undifi_algorithm_steps (b). We call "phantom" those grid-points of the background mesh (shown using dashed circles in Fig. undifi_algorithm_steps (b)) that have been temporarily removed. All cells having at least one phantom node among their vertices are also removed from the background triangulation; these are the cells shown using dashed edges in Fig. undifi_algorithm_steps (b). Further details concerning the criteria used to identify and remove the phantom nodes can be found in [4], [15].

## Local Re-Meshing Around the Shock Front

Following the cell removal step, the background triangulation has been split into two or more disjoint sub-domains, as shown in Fig. undifi_algorithm_steps (c). The hole dug by the fitted front is then re-meshed using a Constrained Delaunay Tessellation (CDT): the edges (triangles in 3D) that make up the fitted discontinuity and the boundary of the hole are both constrained to be part of the final tessellation; this is illustrated in Fig. undifi_algorithm_steps (d). Observe that re-meshing is localized around the discontinuities and, therefore, does not overload the algorithm in terms of CPU cost. Upon completion of this stage, the computational domain is discretized using what we call the "shock-fitting" mesh, which differs from the background mesh only in the neighborhood of the fitted discontinuities. Further details concerning the software used to construct the CDT were given in Sect. The meshing software.

## Calculation of the Unit Vectors Normal to the Shock Front

In order to apply the jump relations, normal (n) and tangent (tau) unit vectors are needed within each pair of grid-points located along the discontinuities, see Fig. fig:algorithm.e. These unit vectors are computed using finite-difference (FD) formulae which involve the coordinates of the shock-point itself and those of its neighboring shock-points. Depending on the local flow regime, it may be necessary to use upwind-biased formulae to avoid geometrical instabilities along the fitted discontinuity. Full details describing how to compute the normals to the discontinuity can be found in [4], [19] for the 2D case and in [4] for the 3D case.

## Solution Update Using the Shock-Capturing Code

Using the shock-fitting mesh as input, a single time step calculation is performed using an unstructured, vertex-centered shock-capturing solver which returns updated nodal values at time t+dt. Since the discontinuities are seen by the shock-capturing code as internal boundaries (of zero thickness) moving with the velocity of the discontinuity, there is no need to modify the spatial discretization scheme already implemented in the PDEs solver to account for the presence of the fitted discontinuities. In practice, the shock-capturing solver is used as a black-box: it receives as input the shock-fitting grid, the nodal values of the solution and grid velocity at time t and returns the updated solution at time t+dt. The solution returned by the shock-capturing solver at time t+dt is however missing some boundary conditions on one or both sides of each discontinuity, depending on whether it is a shock or a contact. These missing pieces of information will be determined in the next step.

## Enforcement of the Jump Relations

The missing pieces of information that are needed to correctly update the solution within all pairs of grid-points located on the discontinuities are obtained by enforcing the R-H jump relations; this also provides the local velocity of the discontinuity along its normal. The R-H jump relations are a set of non-linear algebraic equations that can be solved within each pair of grid-points located along the discontinuities by means of Newton-Raphson's algorithm. In order to match the number of unknowns with the available equations, one or more additional pieces of information are required within both or either of the two sides of the fitted discontinuity, depending on whether this is a shock or a contact discontinuity. These additional pieces of information are obtained from the characteristic formulation of the Euler equations and correspond to those characteristic quantities that are convected towards the discontinuity from the sub-domain that is attached to that side of the discontinuity. Using an upwind-biased discretization within the shock-capturing solver, one can reasonably assume that the spatial and temporal evolution of these characteristic quantities has been correctly computed. Full algorithmic details concerning the practical implementation of the jump relations for shocks and contact discontinuities are reported elsewhere and will not be repeated here: see [9], [15], [16] for the 2D case and [4] for the 3D case. An ad-hoc treatment is required within those special points where different discontinuities interact; this is the case of triple and quadruple points where an impinging shock is reflected from a solid surface, etc. The algorithmic details are described in [9], [16] for the 2D case, whereas the interaction among fitted discontinuities has not yet been dealt with in 3D. This specific issue will be further addressed in Sect. 5.2.

## Shock Displacement

The enforcement of the jump relations provides the speed (w) at which each pair of grid-points located on the discontinuity move along its local normal vector, $n $. The position of the discontinuity at time ... is computed in a Lagrangian manner by displacing all its grid-points, as shown in Fig. fig:algorithm.f where the dashed and solid lines represent the discontinuity at time $t$, resp. ... When simulating steady flows, this can be accomplished using the following first-order-accurate (in time) integration formula: ... which returns the spatial coordinates of the i-th grid-point at time ... The low temporal accuracy of Eq. ? does not affect the spatial accuracy of the steady state solution which only depends on the spatial accuracy of the gas-dynamics solver and that of the tangent and normal unit vectors.

On the contrary, when dealing with unsteady flows, the temporal accuracy of the shock motion has to be the same as that of the spatial discretization, i.e. second order accurate in our case. This can be accomplished using a predictor-corrector type temporal integration scheme, or a Runge-Kutta multi-step scheme. More specifically, the predictor step estimates the position of the shock at time level ... using the explicit Euler scheme:

The shock speed  the normal unit vector ... at time level ... are then computed using the intermediate shock position ... and, finally, the position of each shock point is updated at time level n+1 in the corrector step: ... Figure undifi_algorithm_steps (f) shows that even when the background mesh is fixed in space, the triangular cells that abut on the discontinuity have one of their edges that moves with the discontinuity, thus deforming the cell. This implies that the shock capturing solver used in Step Solution Update Using the Shock-Capturing Code must be capable of handling moving meshes, i.e. it must be capable of solving the governing PDEs written using an Arbitrary Eulerian Lagrangian (ALE) formulation. Finally, the time step ... to be used in Eq. ? to move the shock is chosen in such a way that during the time interval ... the shock will remain within the hole that it has dug in the background mesh. By doing so, none of the grid-points of the shock-fitting mesh will be overcome by the moving discontinuity, as shown in Fig. undifi_algorithm_steps (f).

## Interpolation of the Phantom Nodes

Upon completion of the previous steps, all nodes of the shock-fitting mesh have been updated at time .... The shock-fitting mesh is made up of all the grid-points belonging to the fitted discontinuities and all nodes of the background mesh, except those that have been declared "phantom". Therefore, the nodal values within the phantom nodes have not been updated to time .... However, during the current time step, the discontinuity might have moved sufficiently far away from its previous position, that some of the phantom nodes may re-appear in the shock-fitting mesh at the next time step. It follows that also the nodal values within the phantom nodes need to be updated to time ... This is easily accomplished by transferring the available solution at time ... from the current shock-fitting mesh to the grid-points of the background one, using linear interpolation. Once the phantom nodes have been updated, the shock-fitting mesh used in the current time interval has completed its task and can be removed. At this stage the numerical solution has correctly been updated at time $t+ t$ within all grid-points of the background tessellation and within all pairs of grid-points belonging to the fitted discontinuities. The next time interval can be computed re-starting from the first Step Cell Removal Around the Shock Front of the algorithm.
