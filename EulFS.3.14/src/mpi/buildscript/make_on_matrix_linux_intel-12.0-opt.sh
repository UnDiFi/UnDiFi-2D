module purge
module load caspur
module show compilers/intel
module load compilers/intel
module show libs/acml/4.4.0
module load libs/acml/4.4.0
export PETSC_DIR=/work/ahn/abonfigl/petsc-3.2-p7
export PETSC_ARCH=linux_intel-12.0-opt
export FSPL_DIR=/work/ahn/abonfigl/EulFS.3.2.3

echo Using gcc: `which gcc`
echo
echo Using F95: `which gfortran`
echo

echo PETSC_DIR is: $PETSC_DIR
echo
echo PETSC_ARCH is: $PETSC_ARCH
echo
echo LD_LIBRARY_PATH is: $LD_LIBRARY_PATH
echo

sleep 2

cd $FSPL_DIR/src/mpi

make install

echo
echo "Ricorda di spostare l'eseguibile !!!!!!!!!!!!"
echo
