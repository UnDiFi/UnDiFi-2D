#
# A sh script to compare the sequential and parallel of the code
# differences are sent to ${HOME}/differences
# this script runs mpiseq.sh
#
#
export ITMAX=8
export SNES=newton
export SEQ_CODE=eulfs3.14.0-linux_gnu
export SEQ_CODE=eulfs11.9
export MPI_CODE=${HOME}/bin/linux/peulfs3.14.0-linux_gnu
export MPI_CODE=${HOME}/bin/linux/peulfs11.9
export MPIHOME=${HOME}/mpich
export MPIRUN=${MPIHOME}/bin/mpirun
export MPIRUN=`which mpirun`
#
export NUM_PROCS=2
export DATA_DIR=${HOME}/grids/2D/aerofoils/naca0012/euler/NP2355/
#
# start executing........
#
echo running on `uname -a` `date` > ${HOME}/differences
echo >> ${HOME}/differences
echo "SEQ code is " $SEQ_CODE >> ${HOME}/differences
echo "MPI code is " $MPI_CODE >> ${HOME}/differences
echo >> ${HOME}/differences
#
export TCASE=${HOME}/testcases/2D/compressible.euler/1aA/test/
sh mpiseq.sh ${TCASE}
#
export TCASE=${HOME}/testcases/2D/compressible.euler/1aB/test/
sh mpiseq.sh ${TCASE}
#
# low Mach number case which requires re-starting; re-starting is not readily available in parallel
#
####export TCASE=${HOME}/testcases/2D/compressible.euler/1bB/LDA/test/
####sh mpiseq.sh ${TCASE}
#
export TCASE=${HOME}/testcases/2D/incompressible.euler/1bE/test/
sh mpiseq.sh ${TCASE}
#
# broken for some reason: 11.9 and 0.15.2 work
#
#export TCASE=${HOME}/testcases/2D/incompressible.euler/1bF/test/
#sh mpiseq.sh ${TCASE}
#
export DATA_DIR=${HOME}/grids/2D/aerofoils/naca0012/navier/jcc/
export NUM_PROCS=2
#
export TCASE=${HOME}/testcases/2D/compressible.ns/2aC/jcc/test/
sh mpiseq.sh ${TCASE}
#
export TCASE=${HOME}/testcases/2D/incompressible.ns/2aG/test/
sh mpiseq.sh ${TCASE}
#
export DATA_DIR=${HOME}/grids/3D/wings/sinus/6/
export NUM_PROCS=2
#
export TCASE=${HOME}/testcases/3D/compressible.euler/sinus/6/A/test
sh mpiseq.sh ${TCASE}
#
#
export TCASE=${HOME}/testcases/3D/compressible.euler/sinus/6/B/test
sh mpiseq.sh ${TCASE}
#
# this is used to test the subsonic inflow/outflow b.c.
#
export DATA_DIR=${HOME}/grids/2D/blades/naca/
export NUM_PROCS=2
export TCASE=${HOME}/testcases/2D/compressible.euler/nacaper/test/
sh mpiseq.sh ${TCASE}
#
exit 1
#
# this is used to test scalar problems (advection-diffusion)
#
export DATA_DIR=${HOME}/grids/2D/scalar/advdiff/21x21/
export NUM_PROCS=4
#
export TCASE=${HOME}/testcases/2D/scalar/tcase3
sh mpiseq.sh ${TCASE}
#
# this is used to test scalar problems (advection)
#
export DATA_DIR=${HOME}/grids/2D/scalar/advdiff/41x41/
export NUM_PROCS=4
#
export TCASE=${HOME}/testcases/2D/scalar/tcase5
sh mpiseq.sh ${TCASE}
#
exit 0
#
# this is used to test the turbulent/compressible flows
#
export DATA_DIR=${HOME}/grids/2D/aerofoils/rae2822/navier/jcc/
export NUM_PROCS=2
#
export TCASE=${HOME}/testcases/2D/compressible.ns/rae/
sh mpiseq.sh ${TCASE}
#
# Uh! Oh! Cannot find the partitioned meshes
#
export DATA_DIR=${HOME}/grids/2D/fplate/laminar/Re_5E+5/coarse/
export NUM_PROCS=4
#
export TCASE=${HOME}/testcases/2D/compressible.ns/lbl5e5/
sh mpiseq.sh ${TCASE}
