for dir in `ls -d [QCNMRS]*`

do
   echo $dir
   cd $dir
   ./clean.sh
   pwd
   cd ..
   pwd
done

