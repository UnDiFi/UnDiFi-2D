for dir in seq euler navier-stokes util geometry schemes scalar turbo chemistry mpi
do
	cd $dir
        if [ $dir = "seq" ]
	then
		echo `pwd`
		co -r$1 *.[f,F,c] *.f90 makefile make_intel* make_linux_gnu64*
        elif [ $dir = "mpi" ]
	then
		echo `pwd`
		co -r$1 makefile
        else
		co -r$1 *.[f,F,c] *.f90 Makefile
		echo `pwd`
	fi
	cd ..
done
echo `pwd`
cd ../include
co -r$1 *.com *.h *inc
cd -
