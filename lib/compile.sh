cd libgeometry
./geometry.sh
./geometryf90.sh
cd ..
cd libtriangulation
./triangulation.sh
./triangulationf90.sh
cd ..
cd libtirpc-1.3.1
autoreconf -f -i
./configure
cd src
make 
if [ $? -ne 0 ]
then
        echo "Compiling in libtirpc... failed (make)"
        exit  1
fi
cd .libs
cp libtirpc.a ../../../linux_intel_x86_64/.
cd ../../../../
