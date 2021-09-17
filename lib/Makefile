DIRS=libport \
     libmylib \
     libsparse-blas \
     SPARSKIT2
DIRS2=libfxdr_2.1 \

all:
	./compile.sh
	@for dir in $(DIRS2); \
	do \
	(cd libfxdr_2.1; ./configure; make F77CMD=gfortran; cp *.a ../linux_intel_x86_64; cd .. ); \
	done; \

	@for dir in $(DIRS); \
	do \
	echo Building in $$dir; \
	(cd $$dir; make install; cp *.a ../linux_intel_x86_64; cd .. ); \
	done
