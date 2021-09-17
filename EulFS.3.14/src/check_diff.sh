for dir in seq euler navier-stokes util geometry schemes scalar turbo chemistry mpi
do
	cd $dir
	echo
	echo Current directory is `pwd`
	echo
	if [ $dir = "schemes" ]
	then
 		rcsdiff  *.[f,F,c] Makefile
	elif [ $dir = "chemistry" ]
	then
	 	rcsdiff  *.[f,F] *.f90 Makefile
	elif [ $dir = "mpi" ]
	then
	 	rcsdiff  makefile
	elif [ $dir = "seq" ]
	then
	 	rcsdiff  *.[f,F] makefile
	else
 		rcsdiff  *.[f,F] Makefile
	fi
	cd ..
done
cd ../include
echo
echo Current directory is `pwd` 
echo
rcsdiff  *.com *.h *.inc
cd -
