for dir in `ls -d [QCNMRS]*`

do
   echo $dir
   cd $dir
   # TODO: check parallel conflicts
   ./run.sh -s eulfs -m fitting -f steady | tee run_all_eulfs_fitting_steady_x86.log &
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

