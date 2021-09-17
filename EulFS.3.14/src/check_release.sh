for dir in seq euler navier-stokes util geometry schemes scalar turbo chemistry mpi
do
	cd $dir
	echo
	echo Current directory is `pwd`
	echo
	if [ $dir = "schemes" ]
	then
 		rlog -h  *.[f,F,c] Makefile
	elif [ $dir = "chemistry" ]
	then
	 	rlog -h  *.f *.f90 *.F Makefile
	elif [ $dir = "mpi" ]
	then
	 	rlog -h  makefile
	elif [ $dir = "seq" ]
	then
	 	rlog -h  *.[f,F] makefile
	else
 		rlog -h  *.[f,F] Makefile
	fi
	cd ..
done
cd ../include
echo
echo Current directory is `pwd`
echo
rlog -h  *.com *.h *.inc
cd -
