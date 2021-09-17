FFLAGS	      = -I../../include

include $(PETSC_DIR)/lib/petsc/conf/variables
include $(PETSC_DIR)/lib/petsc/conf/rules

DEST	      = $(FSPL_DIR)/lib/$(PETSC_ARCH)

EXTHDRS	      = ../../include/paramt.h \
		../../include/plasma.h \
		../../include/chem.h \
		../../include/paramchem.h \
		../../include/commonv.inc \
		../../include/commonchem.inc
		
HDRS	      =

INSTALL	      = cp

LIBRARY	      = libchem.a

MAKEFILE      = Makefile

OBJS	      = chemsource.o \
		currentflux.o \
		electcond.o \
		massaction.o \
		ohmsource.o \
		rate.o \
                readArgon.o \
		massactionSU2.o \
		boltzeq4Ar.o \
		iondegree4Ar.o \
		getsigma.o \
		plasflow_inf.o
		
PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = chemsource.f90 \
		currentflux.F \
		electcond.f90 \
		massaction.f90 \
		ohmsource.f90 \
		rate.f90 \
		readArgon.f90 \
		massactionSU2.f90 \
		boltzeq4Ar.f \
		iondegree4Ar.f \
		getsigma.f90 \
		plasflow_inf.f
		
SYSHDRS	      =

all:		$(LIBRARY)

$(LIBRARY):	$(OBJS)
		@echo "Loading $(LIBRARY) ..."
		@ar cru $(LIBRARY) $(OBJS)
		@echo "done"

clobber:;	@rm -f $(OBJS) $(LIBRARY) core tags

depend:;	@mkmf -f $(MAKEFILE)

echo:;		@echo $(HDRS) $(SRCS)

extract:;	@ar x $(DEST)/$(LIBRARY)

index:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(LIBRARY)
		@echo Installing $(LIBRARY) in $(DEST)
		@if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(LIBRARY); $(INSTALL) $(LIBRARY) $(DEST)); fi

print:;		@$(PRINT) $(HDRS) $(SRCS)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)

####include $(PETSC_DIR)/bmake/$(PETSC_ARCH)/base

###
chemsource.o: ../../include/paramt.h ../../include/plasma.h ../../include/constants.h \
        ../../include/paramchem.h ../../include/commonv.inc ../../include/chem.h \
	../../include/commonchem.inc        
massaction.o: ../../include/paramt.h ../../include/paramchem.h \
        ../../include/commonv.inc ../../include/chem.h \
	../../include/commonchem.inc
rate.o: ../../include/paramchem.h ../../include/commonv.inc \
	../../include/commonchem.inc
massactionSU2.o: ../../include/plasma.h ../../include/paramchem.h \
          
