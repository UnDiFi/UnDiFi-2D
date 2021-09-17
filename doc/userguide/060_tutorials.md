\hypertarget{tutorials}{}

# Tutorials \label{chap:tutorials}
This chapter will give a detailed overview of how to set up and run
simulations with UnDiFi-2D. In particular, one steady and one unsteady
testcase will be described. It assumes that some familiar with how to
set the compiler options and how to compile the code. The paths to the
executables are omitted. It is assumed that you either added aliases
for all the executables, or the binary directories to your `$PATH`
variable as described in \ref{sec:installation_directory}.

The provided examples are contained in the `tests` directory. There
are two types of examples:

 * 2D Steady testcases:
   * CircularCylinder-1
   * CoaSHCK
   * MachReflection-1
   * MachReflection-2
   * NACA0012_M080_A0
   * NACA0012_M095_A0_FISHTAIL-1
   * Q1D
   * RegularReflection-1
   * RegularReflection-2
   * SSInteractions1-2
   * SSInteractions2-1
   * SSInteractions2-2
  * 2D Unsteady testcases:
    * ShockVortex

Each testcase can run in shock-capturing and shock-fitting mode. Once
all the simulation options and files are set as explained in Chapter
\ref{chap:workflow}, the simulation can be run by simply launching the
script file **start_run_x86.sh** available in each testcase directory.
The starting step, number of iterations, the gasdynamics solver to
be coupled with UnDiFi-2D, the type of simulation (steady/unsteady)
should be defined. Optionally, output can be redirected to a logfile.

~~~~~~
#                          nbegin, nsteps, eulfs, steady, testcase,         logfile
../../bin/UnDiFi-2D-2D_x86_64 0       501     true   true    "testcase_name" | tee run.log
~~~~~~

**Note:** in future releases of the code, the start_run_x86.sh may
change.

The script **clean.sh** which removes all files and directories is
also available. In the upper level directory there are two scripts to
recursively run and clean all testcases, respectively:

~~~~~~~
./run_all_x86.sh
./clean_all_dir.sh
~~~~~~~

## CircularCylinder-1
This testcase represents a typical blunt body problem which deals
with the high speed flow past a circular cylinder at free-stream Mach
number, $M_1 = 20$. The computational domain surrounds a half circular
cylinder having radius $R = 1$ (see Fig. \ref{blunt1}).

### Creation and visualization of the mesh
The mesh is created using Triangle [@shewchuk1996triangle] mesh
generator, and it is defined by the files na00.1.node, na00.1.poly,
na00.1.edge, na00.1.ele, na00.1.neigh. The files are already present in
the testcase folder but in general the user shall create them according
to the problem's domain.

To visualize the mesh, the following command can be run:

    showme na00.1

and by selecting the appropriate tag it will show the corresponding
mesh entity.

![Computational domain for the CircularCylinder-1 testcase.\label{blunt1}](./images/bluntbody1.png?raw=true)

### Creation of shock file
The solution computed by the unstructured code in shock-capturing mode
is used to initialize the flow-field and to determine the initial
position of the shock front. In the shock-fitting simulation, the
initial upstream state in the shock points is set equal to the free
stream conditions, while the initial downstream state is computed from
the upstream state and the shock slope assuming zero shock speed (
$w=0$ ). The upstream and downstream shock states should be defined
in the input shock file, *sh00.dat* according to a specific format as
explained in the *sh00.dat* file. These values should be in agreement
with those defined in the mesh files otherwise an error may be raised
in the first iterations because the jump conditions which relate
upstream and downstream states are violated.

### Definition of boundary conditions
The next step consists in defining the boundary conditions for the
testcase. Depending of the shock-capturing gasdynamic solver which
will be used, different files are necessary, as described in Chapter
\ref{chap:workflow}.

### Setting input options
The UnDiFi-2D input options are described in Chapter \ref{chap:workflow}.

### Run the simulation
Once completed all the previous steps, it is possible to run the
simulation both in shock-capturing and shock-fitting mode.

#### How to run in shock-capturing mode?
In order to run with the NEO solver, both in shock-capturing and
shock-fitting mode it is necessary to set its basic parameters.

First, move to `NEO_data` and be sure that the `output` is empty.
Then, move to `NEO_data/textinput` and set the main parameters in the
*inputfile-exp.txt* file, as explained in Chapter \ref{chap:workflow}.

To run the simulation in shock-capturing mode, it is sufficient to
modify the settings in `NEO_data/textinput/inputfile-exp.txt` as
desired and directly run the shock-capturing solver from the testcase
directory, by typing:

     ../../bin/CRD_euler

Note that the number of iterations is defined in the field "Maximum
number of time steps" in the `inputfile-exp.txt` file.

To run the simulation in shock-capturing mode with EulFS, launch:

     ../../bin/EulFS_x86_64

#### How to run in shock-fitting mode?
To run the simulation in shock-fitting mode, execute the script:

    ./start_run_x86.sh

The flow-field and the shock position are integrated in time until
steady state is reached. Fig. \ref{blunt1} shows the shock displacement
which occurs between the initial position and the one reached at
steady state. Fig. \ref{blunt2} displays the background mesh and
the computational mesh used during the last iteration. The box on
the right, which shows an enlargement of the near shock region,
allows comparing the two meshes. It can be observed that these are
superimposed everywhere except in the region adjacent to the shock,
where the differences between the background (dashed lines) and the
computational mesh (continuous lines) are due to the addition of
the shock points and edges. Fig. \ref{blunt2} clearly shows that
the re-meshing technique does not significantly increase the number
of points and cells with respect to the background mesh, since the
addition of the shock nodes in the computational mesh is partly
balanced by the removal of the phantom nodes.

![Background and computational mesh. \label{blunt2}](./images/bluntbody2.png?raw=true)

### Visualization of the solution
Solutions are automatically saved in folders named `step#####` each
`IBAK` number iterations, where `IBAK` is defined in the *input.dat*
file. In particular, each folder contains:

 * mesh files of Triangle
 * na99.dat file with the shock nodes
 * sh99.dat file with the updated shock points position and states
 * shocknor.dat file with information about shock normals
 * node velocity file, vel.dat

Depending on the gasdynamic solver adopted, the solution will be
written in different files.

If EulFS is used, in each folder will be present the binary files:

 * file001.dat
 * file002.dat
 * file003.dat

In order to visualize the solution produced by EulFS it is necessary
to post-process the aforementioned files. For this purpose, the scripts
*conv.sh* and *plot0.sh* are provided.

If NEO is used, there will be:

 * NEO mesh file, neogrid.grd
 * NEO initial solution file, vvvv_input.dat
 * NEO solution file, vvvv.dat

Figure \ref{blunt3} shows the comparison between the solutions computed
by the EulFS code working in shock-capturing and shock-fitting mode
on a coarse mesh. Specifically, the coarse grid solution computed in
shock capturing mode is characterized by a very large shock thickness
and by strong spurious disturbances affecting the shock layer region.
These disturbances originate from the shock and are caused by the
mis-alignment between the mesh and the shock. On the contrary, the
shock-fitting solution shows a very neat field. The shock-capturing
solution shows the presence of spurious oscillation even if the shock
thickness has been halved. The fitted solution shows a significant gain
in terms of solution quality with respect the shock-capturing solution.
Of course, the resolution of the shock-capturing solutions could be
significantly improved if local grid-refinement was used in the shock
region, but this typically involves a significant increase in the
number of cells and nodes.

![Comparison of shock-captured and shock-fitted solutions. \label{blunt3}](./images/bluntbody3.png?raw=true)

A quantitative analysis of the shock-fitting solutions was carried
out by comparing the estimates of the shock position and the pressure
distribution computed with the reference solution computed by
Lyubimov and Rusanov [@lyubimov1973gas]. This comparison shows that
the solutions computed by the proposed shock-fitting methodology
are not only in good agreement with the reference solution, but
nearly superimposed. Moreover, this analysis clearly proves that the
coarse grid solution is in practice grid-independent in spite of
the extremely limited number of cells enclosed in the shock layer.
A comparison between the shock-fitting and the shock-capturing
solutions is displayed in Fig. \ref{blunt4} where the normalized
pressure distributions along the line at 45° is plotted. The
shock-capturing solutions are characterized by a finite shock thickness
that significantly affects the distribution even on the fine grid.
Moreover, visible differences between the shock-capturing solutions and
the reference one are also visible in regions far from the shock and
close to the body. On the contrary, the differences between the shock
fitting solutions and the reference one are very small and the solution
is nearly superimposed to the reference one.

![Pressure jump obtained with the captured and fitted solutions. \label{blunt4}](./images/bluntbody6.png?raw=true)

## Shock-Vortex interaction
This test case considers a weak shock-vortex interaction.

The interaction between a shock and a vortex has been frequently
reported in the literature as a tool for understanding the mechanisms
of noise generation due to the interaction between a shock-wave and a
turbulent flow [@grasso2000and].

Figure \ref{SV1} shows the computational domain with the boundary and
the initial conditions. In particular, at the time $t=0$ the field is
characterized by the presence of a normal shock that divides the a
supersonic region from a subsonic region.

![Computational domain for the ShockVortex testcase.\label{SV1}](./images/SVdomain.png?raw=true)

The subsonic field is initialized with the uniform field that is
computed from the steady R-H jump relations for a normal shock wave
with a $M_s=2$ stream in the upstream region. The initial position of
the shock and of the vortex is reported in Fig. \ref{SV1}. Since the
solution for this test is non-stationary and its execution requires
the activation of unsteady mode of the UnDiFi-2D code, this test was
included to verifies all parts and subroutines of the code that ensure
the time accurate computation.

The shock-capturing and shock-fitting simulations were performed
using the NEO solver. The computational domain has been discretized
using a Delaunay triangulation generated with Triangle. The presented
results have been computed on a grid made up of 433664 elements, 217569
points with a $\triangle x = 0.00375$. A qualitative view of the
solutions obtained with the two approaches is given in Fig. \ref{SV2}
and Fig. \ref{SV3}. The pictures show the total enthalpy contours in
the solutions obtained with shock capturing and fitting. In particular,
besides the oscillations related to the approximation of the shock, we
can see clearly that the contours downstream of the discontinuity are
much less smooth in the captured solutions. The fitted computations,
on the other hand, show very nice and smooth contours Further details
about this test case and further simulations concerning the interaction
between shock and vortex can be found in [@campoli2017and].

![Shock-captured solution at $t=0.3$. \label{SV2}](./images/SVct03.png?raw=true)

![Shock-fitted solution at $t=0.3$. \label{SV3}](./images/SVft03.png?raw=true)

All the steps to obtain these solutions are the same as those described
in the previous testcase, thus, they will not be repeated. There is
only one important difference respect to the steady testacases, that
is, for the unsteady testacases, the *timesteps.dat* file must be
present in the directory. This file contains two lines representing the
dt to be used in the predictor-corrector time-accurate integration.

~~~~~~~
 0.0008 # predictor step
 0.0016 # corrector step
~~~~~~~
