#
# A sh script to run the sequential (or parallel) versions
# of the code on different architectures;
# this script runs seq.sh
#
#
export ITMAX=30
export SNES=picard
export SEQ_CODE=eulfs
export MPI_CODE=${HOME}/bin/${PETSC_ARCH}/peulfs
export MPIHOME=/opt/mpich
export MPIRUN=${MPIHOME}/bin/mpirun
#
export NUM_PROCS=2
export DATA_DIR=${HOME}/grids/2D/aerofoils/naca0012/euler/NP2355/
#
# start executing........
#
echo running on `uname -a` `date`
#
#
#
# these are used to test the subsonic, cnst. pressure b.c. 
#
export DATA_DIR=${HOME}/grids/2D/fplate/laminar/Re_5E+5/coarse/
export NUM_PROCS=4
#
export TCASE=${HOME}/testcases/2D/compressible.ns/lbl5e5/
sh seq.sh ${TCASE}
#
export TCASE=${HOME}/testcases/2D/incompressible.ns/lbl5e5/
sh seq.sh ${TCASE}
#
# this is used to test the subsonic inflow/outflow b.c.
#
export DATA_DIR=${HOME}/grids/2D/channels/stanitz/
export TCASE=${HOME}/testcases/2D/incompressible.euler/stanitz/Mach0.69
sh seq.sh ${TCASE}
#
export DATA_DIR=${HOME}/grids/3D/channels/stanitz.euler/
export TCASE=${HOME}/testcases/3D/incompressible.euler/stanitz/
sh seq.sh ${TCASE}
#
