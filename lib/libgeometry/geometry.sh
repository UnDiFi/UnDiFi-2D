#!/bin/bash
#
mkdir temp
cd temp
rm *
./../../../bin/$ARCH/f77split ../geometry.f
#
for FILE in `ls -1 *.f`;
do
  echo compiling $FILE
  gfortran -c -g $FILE >& compiler.txt
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
  rm compiler.txt
done
rm *.f
#
ar qc libgeometry.a *.o
rm *.o
#
mv libgeometry.a ../../linux_intel_x86_64/.

cd ..
rmdir temp
#
echo "Library installed as ../linux_intel_x86_64/libgeometry.a."
