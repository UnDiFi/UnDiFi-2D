#!/bin/bash
#
#ARCH=
mkdir temp
cd temp
rm *
./../../../bin/$ARCH/f90split ../triangulation.f90
#
for FILE in `ls -1 *.f90`;
do
  gfortran -c -g $FILE >& compiler.txt
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
  rm compiler.txt
done
rm *.f90
#
ar qc libtriangulation.a *.o
rm *.o
#
#cp libtriangulation.a ../../linux_intel_x86_64/.
mv -v libtriangulation.a ../../linux_intel_x86_64/libtriangulationf90.a
cd ..
rmdir temp
#
