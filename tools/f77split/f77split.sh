#!/bin/bash
#
gcc -c f77split.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling f77split.c."
  exit
fi
rm compiler.txt
#
gcc f77split.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading f77split.o."
  exit
fi
#
rm f77split.o
#
chmod u+x a.out
mkdir -p ../../bin/$ARCH
mv a.out ../../bin/$ARCH/f77split
#
echo "Executable installed as ../../bin/$ARCH/f77split"
