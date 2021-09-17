mkdir diffdir
dir=util
for i in `ls *.[f,F]`
do
diff $i $FSPL_DIR/../EulFS.3.2.10/src/$dir/$i > diffdir/$i.diff
done
