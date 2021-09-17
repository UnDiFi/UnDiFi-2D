#
# A sh script to diff the convergence histories
# of the code on different architectures;
# this script runs seq.sh
#
echo A sh script to diff the convergence histories
echo of the code on different architectures
#
export ARCH1=linux
export ARCH2=alpha
export MPI_CODE=${HOME}/bin/${PETSC_ARCH}/peulfs
export MPIHOME=/opt/mpich
export MPIRUN=${MPIHOME}/bin/mpirun
#
export NUM_PROCS=2
export DATA_DIR=${HOME}/grids/2D/aerofoils/naca0012/euler/NP2355/
#
# start executing........
#
echo running on `date` > ${HOME}/differences
echo diffs between ${ARCH1} and ${ARCH2} > ${HOME}/differences
#
export TCASE=${HOME}/testcases/2D/1aA
cd ${TCASE}
echo checking ${TCASE}
echo checking ${TCASE} >> ${HOME}/differences
echo diffs are listed below >> ${HOME}/differences
diff convhst.${ARCH1} convhst.${ARCH2} >> ${HOME}/differences
#
export TCASE=${HOME}/testcases/2D/1aB
cd ${TCASE}
echo checking ${TCASE}
echo checking ${TCASE} >> ${HOME}/differences
echo diffs are listed below >> ${HOME}/differences
diff convhst.${ARCH1} convhst.${ARCH2} >> ${HOME}/differences
#
export TCASE=${HOME}/testcases/2D/1bB/LDA
cd ${TCASE}
echo checking ${TCASE}
echo checking ${TCASE} >> ${HOME}/differences
echo diffs are listed below >> ${HOME}/differences
diff convhst.${ARCH1} convhst.${ARCH2} >> ${HOME}/differences
#
export TCASE=${HOME}/testcases/2D/1bE/LDA
cd ${TCASE}
echo checking ${TCASE}
echo checking ${TCASE} >> ${HOME}/differences
echo diffs are listed below >> ${HOME}/differences
diff convhst.${ARCH1} convhst.${ARCH2} >> ${HOME}/differences
#
export TCASE=${HOME}/testcases/2D/1bF/LDA
cd ${TCASE}
echo checking ${TCASE}
echo checking ${TCASE} >> ${HOME}/differences
echo diffs are listed below >> ${HOME}/differences
diff convhst.${ARCH1} convhst.${ARCH2} >> ${HOME}/differences
#
export TCASE=${HOME}/testcases/2D/2aC/jcc/LDA/newton/
cd ${TCASE}
echo checking ${TCASE}
echo checking ${TCASE} >> ${HOME}/differences
echo diffs are listed below >> ${HOME}/differences
diff convhst.${ARCH1} convhst.${ARCH2} >> ${HOME}/differences
#
export TCASE=${HOME}/testcases/2D/2aG/
cd ${TCASE}
echo checking ${TCASE}
echo checking ${TCASE} >> ${HOME}/differences
echo diffs are listed below >> ${HOME}/differences
diff convhst.${ARCH1} convhst.${ARCH2} >> ${HOME}/differences
#
export TCASE=${HOME}/testcases/3D/compressible.euler/sinus/6/
cd ${TCASE}
echo checking ${TCASE}
echo checking ${TCASE} >> ${HOME}/differences
echo diffs are listed below >> ${HOME}/differences
diff convhst.${ARCH1} convhst.${ARCH2} >> ${HOME}/differences
#
export TCASE=${HOME}/testcases/2D/compressible.ns/lbl5e5/
cd ${TCASE}
echo checking ${TCASE}
echo checking ${TCASE} >> ${HOME}/differences
echo diffs are listed below >> ${HOME}/differences
diff convhst.${ARCH1} convhst.${ARCH2} >> ${HOME}/differences
#
#
export TCASE=${HOME}/testcases/2D/incompressible.ns/lbl5e5/
cd ${TCASE}
echo checking ${TCASE}
echo checking ${TCASE} >> ${HOME}/differences
echo diffs are listed below >> ${HOME}/differences
diff convhst.${ARCH1} convhst.${ARCH2} >> ${HOME}/differences
#
#
export TCASE=${HOME}/testcases/2D/incompressible.euler/stanitz/Mach0.69
cd ${TCASE}
echo checking ${TCASE}
echo checking ${TCASE} >> ${HOME}/differences
echo diffs are listed below >> ${HOME}/differences
diff convhst.${ARCH1} convhst.${ARCH2} >> ${HOME}/differences
#
#
export TCASE=${HOME}/testcases/3D/incompressible.euler/stanitz/
cd ${TCASE}
echo checking ${TCASE}
echo checking ${TCASE} >> ${HOME}/differences
echo diffs are listed below >> ${HOME}/differences
diff convhst.${ARCH1} convhst.${ARCH2} >> ${HOME}/differences
#
#
####export DATA_DIR=${HOME}/grids/2D/fplate/laminar/Re_5E+5/coarse/
#
####export TCASE=${HOME}/testcases/2D/compressible.ns/lbl5e5/
#
# this is used to test the subsonic inflow/outflow b.c.
#
####export DATA_DIR=${HOME}/grids/2D/blades/naca/
#
####export TCASE=${HOME}/testcases/2D/compressible.euler/nacaper/test3/
#
# this is used to test the turbulent/compressible flows
#
#
####export TCASE=${HOME}/testcases/2D/compressible.ns/rae/
#
# this is used to test scalar problems (advection-diffusion)
#
#
#####export TCASE=${HOME}/testcases/2D/scalar/tcase3
#
# this is used to test scalar problems (advection)
#
#
#####export TCASE=${HOME}/testcases/2D/scalar/tcase5
