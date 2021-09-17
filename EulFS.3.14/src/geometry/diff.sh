if [ -e diffdir ]
then
	rm -v diffdir/*
else
mkdir diffdir
fi
dir=geometry
for i in `ls *.[f,F]`
do
diff $i $FSPL_DIR/src/$dir/$i > diffdir/$i.diff
done
