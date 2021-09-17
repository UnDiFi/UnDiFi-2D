for dir in `ls -d [QCNMRS]*`

do
   echo $dir
   cd $dir
   # TODO: add path and/or collect logs in one
   ./run.sh -s eulfs -m capturing -f steady | tee run_all_eulfs_capturing_steady_x86.log &
   pwd
   cd ..
   pwd
done

# TODO: check complete
#for dir in `ls -d [QCNMRS]*`
#
#do
#   cd $dir
#   if [ -e step00501 ]
#   then
#       echo test $dir ok
#   else
#       echo test $dir failed
#   fi
#   cd ..
#done

