
export ARCH=`uname -p`
export HOSTTYPE=`uname -p`

echo $ARCH
echo $HOSTTYPE

#echo Compiling in doc
cd doc
./makedoc.sh
echo Done!
cd ..

#echo Compiling in tools
cd tools
./compile.sh
echo Done!
cd ..

#echo Compiling in lib
cd lib
make
if [ $? -ne 0 ]
then
        echo "Compiling in lib ... failed"
        exit  1
fi
echo Done!
cd ..

echo Compiling PETSC 3.14.6
##### Comment these lines if petsc is already downloaded
wget ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.14.6.tar.gz
tar -xvzf petsc-3.14.6.tar.gz
rm petsc-3.14.6.tar.gz
#####
cd petsc-3.14.6
PETSC_DIR1=${PWD}
export PETSC_DIR=${PETSC_DIR1}
./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-mpich --download-fblaslapack
#Extract petsc_arch from newest folder name
PETSC_ARC=$(ls -td -- */ | head -n 1 | cut -d'/' -f1)
export PETSC_ARCH=${PETSC_ARC}
#test petsc-3.14.6
make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} all
if [ $? -ne 0 ]
then
        echo "Compiling PETSC 3.14.6 ... failed"
        exit  1
fi
echo Done!
cd ..

echo Compiling EulFS 3.14
cd EulFS.3.14
FSPL_DIR1=${PWD}
export FSPL_DIR=${FSPL_DIR1}

mkdir ${FSPL_DIR}/lib/${PETSC_ARCH}
echo ${FSPL_DIR}/lib/$PETSC_ARCH

cd src
make
if [ $? -ne 0 ]
then
        echo "Compiling EulFS 3.14 ... failed"
        exit  1
fi
echo Done!
cd ../..

echo Compiling in source_utils
cd source_utils
./compile.sh
if [ $? -ne 0 ]
then
        echo "Compiling in source_utils ... failed"
        exit  1
fi

echo Done!
cd ..

echo
echo Compiling  NEO
cd NEO/src
make clean install
if [ $? -ne 0 ]
then
        echo "Compiling NEO ... failed"
        exit  1
fi
cd ../..

echo Compiling UnDiFi-2D
cd source
make install
if [ $? -ne 0 ]
then
	echo "Compiling UnDiFi-2D ... failed"
	exit  1
fi
echo Done!
cd ..
exit 0
