rm bndflux*
rm convergenza.dat
rm convhst.*
rm file0*
rm fort*
rm check.dat
rm run.log
rm *.log
cd log
rm *
echo ' ' > empty_file
cd ..
rm *.BAK
rm *.plt
rm na?????.*
rm shocknor.dat
rm slip-free.*
rm -r step?????*
rm tecplot.phy
rm timing.000
rm triangle.log
rm vdt010.dat
rm bndry*
rm residual_norm.dat 

# empty NEO dir
cd NEO_data/
cd input/
rm *
echo ' ' > empty_file
cd ..
cd output/
rm *
echo ' ' > empty_file
cd ../../
