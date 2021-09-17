if [ -d diffdir ]
then
	echo "diffdir already exists"
	rm -v diffdir/*.diff
else
mkdir diffdir
fi
echo $FSPL_DIR
dir=seq
for i in `ls *.[f,F] *.f90`
do
diff $i $FSPL_DIR/src/$dir/$i > diffdir/$i.diff
done
