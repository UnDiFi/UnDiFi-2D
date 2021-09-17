for i in `cat log`
do
	echo $i
	co -l $i
	sed s/'X04CAF'/'R8Mat_Print'/ $i > scratchfile; mv -v scratchfile $i
done
