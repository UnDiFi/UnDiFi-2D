\hypertarget{EulFS_tutorial}{}

# EulFS tutorial \label{chap:eultutorial}

This chapter will provide a guide for the user to set up and run
simulations with EulFS, using both a shock capturing and a shock
fitting mode.

## Shock capturing mode
 * **Step 1: how to convert a Triangle mesh into EulFS input files**

As stated in Chap. \ref{workflow}, EulFS code requires the following binary files in input:

 * file001.dat
 * file002.dat
 * file003.dat

 Containing respectively the informations linked to the nodes, the
 grid connectivity and the flow-field variables. They can be generated
 executing the **triangle2dat-NEW_$(HOSTTYPE)** program, provided
 in the *source_utils* directory: the user has to run this code in
 the test-case folder, where the Triangle mesh files (see Chap.
 \ref{chap:mesh} for details) are available.

~~~~~~~
 This code converts triangle files into EulFS files

 Enter fname
~~~~~~~

 At this point, the user has to type the Triangle filename, without any
 extension: the code runs a check on the boundary faces number, then it
 asks for the presence of periodic surfaces:

~~~~~~~
 Are there any periodic surfaces? y/n
~~~~~~~

 If the mesh is periodic, the user enters the two (different) colours
 of the two periodic surfaces: it is also required to enter x/X if
 the corresponding nodes on the two periodic patches have the same x
 coordinate, otherwise enter y/Y (have the same y coordinate). If no
 error messages occur, the code creates the binary files *file???.dat*
 in the same folder of the Triangle mesh.

 * **Step 2: set the runtime options for EulFS in .petscrc file**

 In order to run a numerical simulation, the user has to include in
 this directory also the .petscrc file, which collects the runtime
 options for that particular test-case: all PETSc specific options
 are described in full details in the Users Manual. Some of the EulFS
 options are mandatory, other are specific to the different sets of
 equations the code can deal with. There are also some other options
 that, though necessary to run properly the code, are assigned default
 values whenever the user does not set them explicitly. In particular,
 using the runtime option -colors it is possible to assign a different
 boundary condition to each color, hence to each group of boundary
 faces. Boundary faces are usually grouped together in the meshfile
 by giving the same colour to faces belonging to the same boundary,
 e.g. body faces will be given a different colour for each different
 solid surface, far field will be given another colour, etc. In order
 to assign a specific boundary conditions it is required to fill the
 position° corresponding to the color of the boundary faces with the
 boundary condition code, as explained in the following table:

~~~~~~~
 Colors   Boundary condition
|  0  |   Unused
|  1  |   Supersonic inlet / Dirichlet
|  2  |   Subsonic outlet (constant static pressure)
|  3  |   Supersonic outlet / Neumann
|  4  |   Inviscid wall
|  5  |   Inviscid far field
|  6  |   Viscous wall
|  8  |   Subsonic inlet
~~~~~~~

   ° *The position must be counted after the first -1 value in the -colors line.*

 The solution is saved every *ibak* iterations ( see option -ibak
 in the .petscrc file) in *file010.dat*, which is the EulFS output
 file. Moreover, EulFS code is able to restart a previous numerical
 simulation by reading the solution from file003.dat or file010.dat:
 in this case, the user has to add the following line in the .petscrc
 file:

~~~~~~~
- restart   file003.dat (or file010.dat)
~~~~~~~

  * **Step 3: how to convert EulFS output into Triangle/Tecplot files**

 In order to update the solution stored in the Triangle mesh files,
 EulFS output must be converted into a new *.node* file, which contains
 the grid nodes and the flow-field variables associated with them,
 running **dat2triangle-NEW_$(HOSTTYPE)** in the test-case directory.

~~~~~~~
 This code converts the EulFS code files into triangle files

 Enter triangle file
~~~~~~~

 The user is asked for providing a name for the output *.node*
 file. Then, the program reads the mesh data stored in
 *file001.dat-file002.dat*, the solution contained in *file010.dat* and
 it creates a *filename.node*, using the name entered by the user.

 Post-processing of the solutions computed by EulFS into a
 format suitable for a visualisation is accomplished through
 **fsplot_$(HOSTTYPE)**. It uses the file *inp*, which should be
 present in the same directory, to select the EulFS solution files to
 be converted in a *file012.dat* readable by Tecplot: further details
 about the *inp* syntax can be found in Chapter 7 of the EulFS user
 manual.

## Shock fitting mode
 UnDiFi-2D code runs at each iteration all the steps described above:
 anyway, in order to run EulFS in a shock-fitting mode, the following
 options must be present in .petscrc file:

 * Only an iteration must be performed:

~~~~~~~
- itmax 1
~~~~~~~

 * EulFS computation must restart from file003.dat, where the solution computed in the previous iteration is stored.

~~~~~~~
- restart   file003.dat
~~~~~~~

* Shock edges, which are associated to color 10 in the computational
mesh, need no particular boundary condition in EulFS: thus, the option
-colors in .petscrc file must contain a 0 in the 10-th position.

~~~~~~~
- colors   -1,.,.,.,.,.,.,.,.,.,0,.,.,.
~~~~~~~
