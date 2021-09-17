\hypertarget{workflow}{}

# Workflow \label{chap:workflow}
In this chapter, the complete process of setting up a simulation with
UnDiFi-2D is detailed.

## What happens when running compile_all.sh? \label{sec:compileroptions}
Let's see what happens when the script **compile_all.sh** is executed.

First, we generate the documentation:

~~~~~~~
echo Compiling in doc
cd doc
./makedoc.sh
echo Done!
cd ..
~~~~~~~

Then, f77split and f90split are compiled. They split the modules of a
fortran file into separate files. A "module" is a blockdata, function,
module, program, or subroutine program subunit.

~~~~~~~
echo Compiling in tools
cd tools
./compile.sh
echo Done!
cd ..
~~~~~~~

Then, several libraries are compiled, specifically:

 * libfxdr_2.1
 * libport
 * libmylib
 * libsparse-blas
 * SPARSKIT2
 * libgeometry
 * libtriangulation
 * libtirpc

~~~~~~~
echo Compiling in lib
cd lib
make
echo Done!
cd ..
~~~~~~~

Then, the [PETSC](https://www.mcs.anl.gov/petsc/index.html) library
is downloaded and compiled. **Note**: it is recommended to use the
proposed version of the library in order to avoid issues with the
compilation of EulFS. It may be updated in future releases.

The C and Fortran compilers in the configuration settings may be
changed. Tested compilers are Intel and GNU pre-installed version of
MPI and LAPACK libraries may be used. In this case, the PATH should be
provided.

~~~~~~~
echo Compiling PETSC 3.14.6
wget ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.14.6.tar.gz
tar -xvzf petsc-3.14.6.tar.gz
rm petsc-3.14.6.tar.gz
#
cd petsc-3.14.6
PETSC_DIR1=${PWD}
export PETSC_DIR=${PETSC_DIR1}
./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-mpich --download-fblaslapack
#Extract petsc_arch from newest folder name
PETSC_ARC=$(ls -td -- */ | head -n 1 | cut -d'/' -f1)
export PETSC_ARCH=${PETSC_ARC}
#test petsc-3.14.6
make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} all
echo Done!
cd ..
~~~~~~~

The gasdynamic solver, EulFS, is then compiled and installed:

~~~~~~~
echo Compiling EulFS 3.14
cd EulFS.3.14

FSPL_DIR1=${PWD}
export FSPL_DIR=${FSPL_DIR1}

mkdir ${FSPL_DIR}/lib/${PETSC_ARCH}
echo ${FSPL_DIR}/lib/$PETSC_ARCH

cd src
make 
echo Done!
cd ../..
~~~~~~~

As explained in the previous section, in addition to the
shock-capturing gasdynamic solvers, EulFS and NEO, there are several
auxiliary codes, contained in the `source_utils` directory, which
are called during the steps of the UnDiFi-2D algorithm or during the
post-processing phase. In particular, the following packages are
compiled and installed:

* dat2triangle: reads files in EulFS format then writes two input files for Triangle, .node and .poly
* dat2paraview: converts both Neo and EulFS output into Paraview files (*.dat)

* Grid_0: creates the neogrid0.grd file
* Na00x2vvvv: writes the file vvvv_input.dat
* Na_creation: creates the triangle mesh files for pre-set testcases
<!-- * NEO_source: -->
* NEO2triangle: reads files in NEO format then writes two input files for Triangle, .node and .poly
* triangle2dat: converts a Triangle mesh to a .dat file format. The input file for Triangle can be created using dat2triangle
* Triangle2grd: converts triangle format to .grd format (NEO)
* Triangle: mesh generator for UnDiFi-2D
<!-- * Tecplot: post-processing output generator for UnDiFi-2D
* fsplot: performs post-processing for EulFS
-->

**Note**: the Tecplot visualizer, may not be necessary as solutions can
be outputted in ASCII format and post-processed by Paraview or Visit,
open-source softwares.

~~~~~~~
echo Compiling in source_utils
cd source_utils
./compile.sh
echo Done!
cd ..
~~~~~~~

Finally, UnDiFi-2D is compiled and installed.

~~~~~~~
echo Compiling in source
cd source
make install
echo Done!
cd ..
~~~~~~~

## Mesh generation/conversion with Triangle
The shock is treated as an internal boundary by the unstructured
gasdynamic solver so that no modifications in the computational kernel
of the CFD code are required. Nonetheless, the shock points can freely
move over the mesh points of the background grid, as it was done in
the floating-shock fitting approaches proposed in the past. To achieve
this, the computational mesh is locally re-meshed at each time-step
using the background grid and the shock edges and it differs from the
background mesh only in the neighbourhood of the shock front, thus
keeping the overall number of grid points to a minimum.

For the present release version of UnDiFi-2D, the preferred mesh generator
is [Triangle](https://www.cs.cmu.edu/~quake/triangle.html), which will
be described in Chapter \ref{chap:mesh}.

## Solver settings input files
Before setting up a simulation, depending on the gasdynamic solver
used, few parameters must be set in the input files.

### NEO
In order to use UnDiFi-2D coupled with NEO, two files must be provided:
**type.dat** and **inputfile-exp.txt**.

The ascii file type.dat contains the definition of boundary
conditions needed by NEO and it is read by the auxiliary utility
converter program Triangle2grd.

~~~~~~~
Face  Type
1      1
2     -1
3      3
4     -1
5     -1
6     -1
7     -1
8     -1
9     -1
10     2

# boundary #1 is the inflow
# boundary #2 is the upper outlow
# boundary #3 is the wall
# boundary #4 is the upper outlow
#
# This file is read by Triangle2grd
#
# bndry code 1 is supersonic
# bndry code -1 is outflow
# bndry code 3 is inviscid wall
~~~~~~~

**Note**: it is implicitely assumed that the computational domain is
made up at most of 10 edges. If more that 10 edges are present, the
source code in `Triangle2grd` should be modified accordingly, and
faces added in this file.

**Note**: in the case that other boundary conditions are necessary,
they should be implemented in NEO, an additional Type flag added in the
source code of `Triangle2grd` and this file, as well.

**Note**: the type.dat file is not required to use NEO in
shock-capturing mode.

The inputfile-exp.txt file defines all the other settings required by
NEO, as detailed below. It must be located in the `NEO_data/textinput`
folder.

~~~~~~~
          NEO's basic parameters

                 PHYSICS

Model                                       0
Initial state                               0
Initial solution file              vvvv_input
Problem                                     0
Gamma                                     1.4
Zero density                         0.000001

                 GEOMETRY

External grid file                    neogrid

                 NUMERICS

0(LDA),1(LLF),2(LDAN),3(LWLF),4(LW),5(N)    5
Lumping type: 0 (selective),1 (global)      0
Shield factor                             0.2

                 ALGEBRA

Convergence treshold                     -3.0
Convergence limit                        -13.
Check Steady State                          0
Maximum number of iterations                2
Convergence variable (0->Neq.s-1)          -1

                  STOP

Final time                           1000000.
Maximum number of time steps                1
Stop                                        0

                 OUTPUT

Movie (0/1)                                 0
Steps to save                             100
Information Frequency                       1
Output format  (0/1)                        1
Output file                              vvvv
~~~~~~~

For unsteady simulations the file timesteps.dat must be provided
for the NEO solver. It contains the time step values used in the
Predictor-Corrector time accurate integration method.

~~~~~~~
0.002 # predictor step dt
0.004 # corrector step dt
~~~~~~~

**Note:** this file is required by the current version of NEO but it
may not be necessary in feature releases due to a change in the time
accurate integration method.

### EulFS
In order to use UnDiFi-2D coupled with EulFS, one file must be provided:
**.petscrc**.

~~~~~~~
-equation Euler
-fluid compressible
-preconditioning No
-nondimensionalisation external

-scalar_scheme NL
-matrix_scheme B

#-scalar_scheme N
#-matrix_scheme N

#-scalar_scheme LDA
#-matrix_scheme LDA

-itmax 1
-ibak 500
-cfl 0.5
-timestepping explicit
-timestep local
-cflmax 1.
-linearization picard

# flow conditions
#
# Mach = 0.63
# 2 degrees of incidence

-freestream_Mach_number 20.0
-flow_angles -1.,0.,0.
-bc_type weak
-colors -1,5,5,4,5,-1,-1,-1,-1,-1,0,-1,-1,-1
-data_dir ./
-restart_file file003.dat
~~~~~~~
The syntax used in this file is defined by the PETSC library.
The reader is referred to its documentation for further details.

All other options are set in the following parameter files.

## UnDiFi-2D setup parameter input files
In addition to the files required by the solvers, two other files are
necessary for UnDiFi-2D: **input.dat** and **sh00.dat**.

The ascii file input.dat contains the main options to set and tune
UnDiFi-2D. A template of this file is as following:

~~~~~~~
0.20d-9   ! EPS
0.30      ! SNDMIN
0.10      ! DXCELL
0.9       ! SHRELAX
50        ! IBAK
1.40d+0   ! GA
0         ! number of additional hole points
0         ! number of periodic boundary
0         ! filter on discontinuity speeds  [0-1.0]  0=disactive
0         ! iteration of mesh topology freezing (must be equal to a backup iteration)
~~~~~~~

 `EPS` represents the minimum relative distance between shock edges \
 `SNDMIN`	\
 `DXCELL`	represents an estimate of the cell edge length \
 `SHRELAX` allows for a relaxation from/to the captured solution \
 `IBAK`	sets the solution saving rate \
 `GA`	is the gas constant

The sh00.dat file contains the structure of the discontinuities present
in the flow field at the beginning of the simulation. An example is
reported below.

~~~~~~~
1      # number of discontinuities
22 'S' # number and type of points forming the discontinuity.
  2  0.30 0.1000E+01 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
  3  0.30 0.9523E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
  4  0.30 0.9047E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
  5  0.30 0.8571E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
  6  0.30 0.8095E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
  7  0.30 0.7619E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
  8  0.30 0.7142E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
  9  0.30 0.6666E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
 10  0.30 0.6190E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
 11  0.30 0.5714E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
 12  0.30 0.5238E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
 13  0.30 0.4761E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
 14  0.30 0.4285E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
 15  0.30 0.3809E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
 16  0.30 0.3333E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
 17  0.30 0.2857E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
 18  0.30 0.2380E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
 19  0.30 0.1904E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
 20  0.30 0.1428E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
 21  0.30 0.9523E-01 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
 22  0.30 0.4761E-01 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
 23  0.30 0.0000E+00 0.1639E+01 0.5944E+01 0.1204E+01 0.0 0.1183E+01 0.2958E+01 0.0 0
 24 2        # number of special points
 25 'WPNRX'  # type of special points: wall point normal ...
 26 1 1      # shock 1, ending point 1
 27 'WPNRX'  #
 28 1 2      # shock 1, ending point 2
~~~~~~~

The data points correspond, respectively, to the X and Y coordinates
(in 2D) and to the downstream and upstream shock states, in term of
Roe variables Z defined as Z = $\sqrt{\rho}$ (1, H, u), where H is the
total enthalpy and u the velocity vector.

## Converters
As stated elsewhere, in order to couple the UnDiFi-2D algorithm code
with gas-dynamic solvers, it is necessary to provide "converters" which
allow data transfer between possibly different formats.

It is important to note that each converter is called during the
normal running mode of the UnDiFi-2D code but it can also be used in
"stand-alone mode", as described below. This capability can be useful
for debugging or testing purposes or when creating or restarting a
simulation. In general, a simple script is provided to compile them in
the corresponding directory.

### Na00x2vvvv
This converter writes the triangle mesh files into the NEO solution
data format file, vvvv.dat.

To compile it, just run, for example:

~~~~~~~
gfortran -o na2vvvv main_n2v.f90
~~~~~~~

By default, it read the na00.?.node file and writes the
vvvv_input.dat file in the `NEO_data/output` directory.

### Na_creation
This program creates na00.node and na00.poly files for a steady planar
shock. It could be easily adjusted to generate other testcases.

To compile it, just run, for example:

~~~~~~~
gfortran -o exe_na physical_values.f90 reading_mod.f90 writing_mod.f90 main.f90
~~~~~~~

### Grid_0
This program creates the neogrid0.grd which is the grid file used by NEO when
running in capturing mode.

To compile it, just run, for example:

~~~~~~~
gfortran -o neogrid0 main.f90
~~~~~~~

### NEO2triangle
This converter reads the NEO grid file neogrid.grd in
`NEO_data/input` and the solution file vvvv.dat at the current time
step, in `NEO_data/output` and writes the corresponding triangle mesh
file na00.?.1.node and na00.?.1.poly file. Note that the .poly file
does not change respect to the previous step as far as the element
connectivity does not change.

To compile it, just run, for example:

~~~~~~~
gfortran -o NEO2triangle main_d2t.f90
~~~~~~~

### Triangle2grd
This program converts the mesh triangle file into the NEO (.grd) mesh
format. It requires the triangle mesh files `.node, .ele, edge, .neigh,
.poly` and `type.dat` and writes the `neogrid.grd`.

To compile it, just run, for example:

~~~~~~~
gfortran -o triangle2grd boundary_edge_mod.f90 r_files_mod.f90 main.f90
~~~~~~~

### dat2triangle
This program reads files in EulFS format and then writes two input
files for triangle: *.node and *.poly. It uses the file *inp* which should
be present in the same directory to select the EulFS solution files,
typically file001.dat, file002.dat, file003.dat.

### triangle2dat
This program converts a triangle mesh to a *.dat* file format used by
EulFS. The input file for triangle can be created using dat2triangle.
This program needs .edge .neigh and .poly files, hence run:

~~~~~~
triangle -n -e -p filename.poly
~~~~~~

where:

* -n outputs (to a .neigh file) a list of triangles neighboring each triangle.
* -e outputs (to an .edge file) a list of edges of the triangulation.
* -p triangulates a Planar Straight Line Graph (.poly file).

For a detailed explanation about all the flags available in Triangle, the user
is referred to its main page.

## Simulation
After the mesh generation, compilation of the binary and setup
of the parameter files, the code can be executed by running the
**start_run_x86.sh** script file present in each testcase directory.

For convenience, a script file is made available, **run.sh** which allows
to launch any simulation. Its usage is the following:

~~~~~~
./run.sh -s Solver -m Mode -f Flow
	-s Solver: neo or eulfs
	-m Mode: capturing or fitting
	-f Flow: steady or unsteady
~~~~~~

Each `IBAK` iterations, where the parameter `IBAK` has been defined in
the input.dat file, the solution is saved in a newly created folder
named `step#####` corresponding to the saving rate.

Depending on the gasdynamic solver coupled with UnDiFi-2D,
different files will be present in the solution folders.

The following files will be always present in the solution folders.

 * na0000?.1.edge:
 * na0000?.1.ele:
 * na0000?.1.neigh:
 * na0000?.1.node:
 * na0000?.1.node.BAK:
 * na0000?.1.poly:
 * na0000?.node:
 * na0000?.poly:
 * na99.node:
 * sh99.dat:
 * shocknor.dat:

If EulFS it is used, the following additional files will be present:

 * file001.dat:
 * file002.dat:
 * file003.dat:
 * file010.dat:

If NEO it is used, the following additional files will be present:

 * neogrid.grd: grid file
 * vvvv.dat: flow field solution
 * vvvv_input.dat: initial flow field solution

## Post-processing \label{sec:postproc}
Solution files generated by NEO are readily available in each
solution step folder, step?????. At the moment, the preferred
visualizer is Tecplot but Paraview/Visit can also be selected in the
inputfile-exp.txt file.

Post-processing of the solutions computed by EulFS into a format
suitable for a visualisation is accomplished through fsplot which
reads run-time options from a control file named *inp* located in the
same directory. More information about fsplot can be found in Chapter 7
of the EulFS user manual present in its `doc` folder.

To make the conversion process easier, a script file, *plot.sh* is made
available in the testcase directory: it requires to specify the solver 
used for the computation (neo or eulfs option), as shown below.
~~~~~~~
 ./plot.sh -s  solver
~~~~~~~
Conversion instructions are briefly reported below: in particular, Paraview users
have to uncomment the lines which refer to dat2paraview converter: in this case, the code 
will automatically detect whether an EulFS or Neo output file is stored in the step????? 
folder. Then a step?????.dat file will be created in the test-case directory and it can be 
read by Paraview using the Tecplot reader option.

~~~~~~~
 for dir in `ls -d step?????`

do
   echo $dir
   cp inp $dir
   cd $dir
   if [ "$Solver" = "eulfs" ]
   then
     ../../../bin/fsplot-$HOSTTYPE
      preplot file012.dat $dir-eulfs.plt
   elif [ "$Solver" = "neo" ]
   then
      preplot vvvv.dat $dir-neo.plt
   fi
   
###### FOR PARAVIEW USERS #######
#      ../../../bin/dat2paraview
#        cp paraview.dat ../$dir.dat  
##############################		
   
   mv $dir*.plt ..
   cd ..
done

~~~~~~~

This conversion step can be directly incorporated into the run.sh script
to make it automatic.
