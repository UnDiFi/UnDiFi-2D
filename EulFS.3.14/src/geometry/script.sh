for i in *.[f,F]
do
	rlog -h $i > log
	ver1=`grep head log | cut -f2 -d:`
	ver2=`grep $1 log | cut -f2 -d:`
	echo "Now checking" $i
	echo "    release" $1 "has version" $ver2 " head is " $ver1
	if [ $ver1 != $ver2 ]
        then
	echo "smthg. wrong"; exit 1
        fi
done
