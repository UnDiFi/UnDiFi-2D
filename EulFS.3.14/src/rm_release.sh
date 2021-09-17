for dir in seq euler navier-stokes util geometry schemes scalar turbo chemistry mpi
do
	cd $dir
	echo `pwd`
	if [ $dir = "schemes" ]
	then
 		rcs -n$1 *.[f,F,c] Makefile
	elif [ $dir = "chemistry" ]
	then
	 	rcs -n$1 *.[f,f90] Makefile
	elif [ $dir = "mpi" ]
	then
	 	rcs -n$1 makefile
	elif [ $dir = "seq" ]
	then
	 	rcs -n$1 *.[f,F] makefile
	else
 		rcs -n$1 *.[f,F] Makefile
	fi
	cd ..
done
cd ../include
echo `pwd`
rcs -n$1 *.com *.h *.inc
cd -
