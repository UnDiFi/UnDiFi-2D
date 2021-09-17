FFLAGS = 
FC = g77

#FC = cf77
#FFLAGS = -Wf"-dp"

FILES = rinfo1.o dinfo13.o

info1.ex: $(FILES)  ../libskit.a 
	$(FC) $(FFLAGS) -o info1.ex $(FILES) ../libskit.a 

clean:
	rm -f *.o *.ex core *.trace

../libskit.a:
	(cd ..; $(MAKE) $(MAKEFLAGS) libskit.a)

