export PETSC_DIR=/usr/local/src/petsc/petsc-2.3.3-p8
arch1=linux_g95-opt
arch2=linux_intel-opt
logfile=$FSPL_DIR/testexamples/$arch1\_vs\_$arch2.log
echo Using logfile: $logfile
rm -f $logfile
echo
#
# outer loop on picard/newton
#
#for snes in picard newton
for snes in newton
do
#
# inner loop on compiler version
#
  for petsc_arch in $arch1 $arch2
  do
     eulfs=eulfs11.13-$petsc_arch
     make testmachines "SNES=$snes" "PETSC_ARCH=$petsc_arch" "SEQ_CODE=$eulfs" "SCALAR_SCHEME=LDA" "MATRIX_SCHEME=LDA" "LOGFILE=$logfile"
  done # end loop on compiler versions
  echo "----------------------------------------------" >> $logfile
  echo "----------------------------------------------" >> $logfile
  echo
  echo "linearization is set to " $snes >> $logfile
  echo
  echo "----------------------------------------------" >> $logfile
  echo "----------------------------------------------" >> $logfile
  make diffmachines "ARCH1=linux_g95-opt" "ARCH2=linux_intel-opt" "LOGFILE=$logfile"
done # end loop on snes
#
# outer loop on schemes
#
echo "**********************************************" >> $logfile
echo "**********************************************" >> $logfile
echo
echo "Now checking differences btw matrix schemes " >> $logfile
echo
echo "**********************************************" >> $logfile
echo "**********************************************" >> $logfile
snes=picard
for matrixscheme in LW LW2 LDA2
do
echo
echo
echo matrix scheme is now set to $matrixscheme
echo
echo
#
# inner loop on compiler version
#
  for petsc_arch in $arch1 $arch2
  do
     eulfs=eulfs11.13-$petsc_arch
     make testmachines "SNES=$snes" "PETSC_ARCH=$petsc_arch" "SEQ_CODE=$eulfs" "SCALAR_SCHEME=LDA" "MATRIX_SCHEME=$matrixscheme" "LOGFILE=$logfile"
  done
  echo "----------------------------------------------" >> $logfile
  echo "----------------------------------------------" >> $logfile
  echo
  echo "matrix scheme is set to " $matrixscheme >> $logfile
  echo
  echo "----------------------------------------------" >> $logfile
  echo "----------------------------------------------" >> $logfile
  make diffmachines "ARCH1=linux_g95-opt" "ARCH2=linux_intel-opt" "LOGFILE=$logfile"
done # end loop on matrix schemes
