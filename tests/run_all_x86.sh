for dir in `ls -d [QCNMRS]*`

do
   echo $dir
   cd $dir
   ./start_run_x86.sh
   pwd
   cd ..
   pwd
done

for dir in `ls -d [QCNMRS]*`

do
   cd $dir 
   if [ -e step00501 ]  
   then 
       echo test $dir ok 
   else
       echo test $dir failed 
   fi
   cd ..
done

