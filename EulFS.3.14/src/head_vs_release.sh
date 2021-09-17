#
# Here we check that the head coicides with the version in the releaseX_Y_Z
#
# ./head_vs_release.sh releaseX_Y_Z
#
#
check_head_vs_release ()
#
# first argument is the file
# second argument is the release
#
#
{
	rm log
	rlog -h $1 >& log
	head=`grep head log|cut -f2 -d:` 
	head_release=`grep $2 log|cut -f2 -d:` 
	echo "Now processing " $1
	if [[ $head != $head_release ]]
	then
		echo "Uh! Oh! the current head" $head "is NOT in the release " $1 $head_release
		echo "check" $1 "in" `pwd`
		return -999
	else
		return 0
	fi
}
for dir in seq euler navier-stokes util geometry schemes scalar turbo chemistry mpi
do
	cd $dir
	echo
	echo `pwd`
	echo
	if [ $dir = "schemes" ]
	then
		for i in *.[f,F,c] Makefile
		do
		check_head_vs_release $i $1
		done
	elif [ $dir = "chemistry" ]
	then
		for i in *.[F,f] *.f90 Makefile
		do
		check_head_vs_release $i $1
		done
	elif [ $dir = "mpi" ]
	then
		for i in makefile
		do
		check_head_vs_release $i $1
		done
	elif [ $dir = "seq" ]
	then
		for i in  *.[f,F] makefile make_*
		do
		check_head_vs_release $i $1
		done
	else
		for i in  *.[f,F] Makefile
		do
		check_head_vs_release $i $1
		done
	fi
	cd ..
done
cd ../include
echo `pwd`
		for i in *.com *.h *.inc
		do
		check_head_vs_release $i $1
		done
