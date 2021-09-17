rm bndflux* bndry0?.dat
rm convergenza.dat
rm convhst.*
rm fort*
rm ale.log
rm check.dat
rm *.log
rm integral.000

rm grid*

rm file0* fold0*

cd log
rm *
cd ..

cd Tec
rm vvvv*
cd ..

cp timesteps.dat.BKP timesteps.dat

cd vel
rm *
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
rm residual_norm.dat

cd NEO_data
cd output
rm vvvv* 
rm *.plt
cd ..

cd textinput
mv inputfile-exp.txt.BAK inputfile-exp.txt
cd ..

cd input
rm vel.dat
rm neogrid0.grd
rm neogrid.grd
