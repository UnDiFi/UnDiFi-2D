for dir in `ls -d step?????`

do
   echo $dir
   cp inp $dir
   cd $dir
   ../../../bin/fsplot-$HOSTTYPE
   preplot file012.dat $dir.plt
   mv $dir.plt ..
   cd ..
done
