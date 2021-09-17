for i in `cat log`
do
	echo $i
	sed s/'X04EAF'/'I4Mat_Print'/ $i > scratchfile; mv -v scratchfile $i
done
