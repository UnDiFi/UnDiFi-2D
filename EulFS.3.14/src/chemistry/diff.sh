mkdir diffdir
echo $FSPL_DIR
dir=chemistry
for i in `ls *.[f] *.f90`
do
diff $i $FSPL_DIR/src/chemistry/$i > diffdir/$i.diff
done
