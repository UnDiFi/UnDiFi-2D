if (test -s $1.node && test -s $1.elem)
then
	echo Processing file : $1
else
	echo Empty or non existing file : $1
	exit
fi
NPOIN=`cat $1.node | wc --lines`
echo There seem to be $NPOIN meshpoints in $1.node
NELEM=`cat $1.elem | wc --lines`
echo There seem to be $NELEM cells in $1.elem
echo TITLE = '"STAR generated mesh"' > $1.dat
echo VARIABLES = '"X","Y","Z"' >> $1.dat
echo ZONE T='" 3"', F=FEPOINT,ET=BRICK      ,N= $NPOIN,E= $NELEM >> $1.dat
cut -c9- $1.node >> $1.dat
cut -c1-70 $1.elem >> $1.dat
echo Tecplot ASCII file $1.dat has been created
preplot $1.dat
echo Tecplot binary file $1.plt has been created
