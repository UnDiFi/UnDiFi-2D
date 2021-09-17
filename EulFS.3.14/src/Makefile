DIRS=	euler\
	geometry\
	chemistry\
	util\
	navier-stokes\
	scalar\
	schemes\
	seq\
	turbo\
	mpi\
	seq
clean: 
	@echo
	@echo FSPL_DIR is currently set to  $(FSPL_DIR)
	@echo
	@echo PETSC_ARCH is currently set to  $(PETSC_ARCH)
	@echo
	@echo PETSC_DIR is currently set to  $(PETSC_DIR)
	@echo
	sleep 1
	@for dir in $(DIRS) ;\
          do \
          echo cleaning $$dir ;\
          (cd $$dir; make -s clean) ;\
          done; rm -f ../lib/$(PETSC_ARCH)/*.a seq/eulfs mpi/peulfs;\
          find . -name "*.a" -exec rm {} \; 
	@for dir in $(DIRS) ;\
          do \
	      (cd $$dir; make install; cd ..) ;\
	  done



