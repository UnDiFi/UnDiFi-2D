# rcsdiff *.[f,F] makefile make_*
# rcsdiff -q --rcs *.[f,F] makefile make_*
rm -f log
rm -f lista
echo `date` > log
echo `date` > lista
for i in *.[f,F] Makefile
do
	rcsdiff $i
	retval=$?
	if [ $retval -eq 1 ]
	then
           echo $i >> log
	   rcsdiff $i >> log 
	   echo "There are difference with" $i $retval
	   echo $i >> lista
	fi
done
