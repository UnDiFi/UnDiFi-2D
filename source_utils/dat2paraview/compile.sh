echo compiling dat2paraview
gfortran -c dat2paraview.f
gfortran -o dat2paraview dat2paraview.o
cp dat2paraview ../../bin/.