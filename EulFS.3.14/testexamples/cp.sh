for i in file010.dat fold010.dat 
do 
	if [ -e $i ] 
	then
		cp $i $i.BAK
	else
		echo 'There is no' $i 'to be copied' ; \
	fi
	done
