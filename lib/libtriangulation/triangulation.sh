#!/bin/bash
#
#ARCH=
mkdir temp
cd temp
./../../../bin/$ARCH/f77split ../triangulation.f
#
for FILE in `ls -1 *.f`;
do
  gfortran -c -g $FILE >& compiler.txt
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
  rm compiler.txt
done
rm *.f
#
ar qc libtriangulation.a *.o
rm *.o
#
cp  libtriangulation.a ../../linux_intel_x86_64/.
cd ..
rmdir temp
#
