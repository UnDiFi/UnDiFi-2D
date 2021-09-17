for i in readat.F weakbc.F solzne.F
do
echo $i
sed s/'include\/finclude'/'finclude'/ $i > tmpfile; mv tmpfile $i
done
