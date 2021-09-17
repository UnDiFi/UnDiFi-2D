if [ -d diffdir ]
then
	echo "diffdir already exists"
	rm -f diffdir/*.diff
else
mkdir diffdir
fi
echo $FSPL_DIR
dir=include
for i in `ls *.h *.com *.inc`
do
diff $i $FSPL_DIR/$dir/$i > diffdir/$i.diff
done
