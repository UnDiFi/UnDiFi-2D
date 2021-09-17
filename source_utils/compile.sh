echo
echo Compiling triangle2dat
cd triangle2dat
#make clean
make clean install
if [ $? -ne 0 ]
then
        echo "Compiling triangle2dat  ... failed"
        exit  1
fi

cd ..

echo Compiling in dir dat2triangle
cd dat2triangle
#make clean
make clean install
if [ $? -ne 0 ]
then
        echo "Compiling dat2triangle  ... failed"
        exit  1
fi


cd ..

cd fsplot/plot
echo "Compiling fsplot "
echo "About to compile fsplot; check the log file, if needed"
echo
sleep 1
#make clean
make clean install
if [ $? -ne 0 ]
then
        echo "Compiling fsplot  ... failed"
        exit  1
fi

cd ../../
echo
#
echo
echo Compiling triangle-v1.6
cd triangle-v1.6
make install
if [ $? -ne 0 ]
then
        echo "Compiling triangle-v1.6  ... failed"
        exit  1
fi


cd ..
#
echo
echo Compiling NEO2triangle
cd NEO2triangle
#./compile.sh
make clean install
if [ $? -ne 0 ]
then
        echo "Compiling NEO2triangle  ... failed"
        exit  1
fi
cd ..
#
echo
echo compiling in Na00x2vvvv
cd Na00x2vvvv
#./compile.sh
make clean install
if [ $? -ne 0 ]
then
        echo "Compiling Na00x2vvvv  ... failed"
        exit  1
fi
cd ..
#
echo
echo Compiling  Triangle2grd
cd Triangle2grd
#./compile.sh
make clean install
if [ $? -ne 0 ]
then
        echo "Compiling Triangle2grd  ... failed"
        exit  1
fi
cd ..
#
echo
echo Compiling neogrid0
cd Grid_0
#./compile.sh
make clean install
if [ $? -ne 0 ]
then
        echo "Compiling neogrid0  ... failed"
        exit  1
fi
cd ..
#
echo
echo Compiling in Na_creation
cd Na_creation
./compile.sh
cd ..

echo Compiling in Dat2paraview
cd dat2paraview
./compile.sh
cd ..
#

