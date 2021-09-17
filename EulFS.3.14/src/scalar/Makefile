include $(PETSC_DIR)/lib/petsc/conf/variables
include $(PETSC_DIR)/lib/petsc/conf/rules

FFLAGS	      = -I../../include

DEST	      = $(FSPL_DIR)/lib/$(PETSC_ARCH)/.

EXTHDRS	      = ../../include/constants \
		../../include/flags.com \
		../../include/io.com \
		../../include/paramt.h \
		../../include/scalar.com \
		../../include/three \
		../../include/visco.com

HDRS	      =

INSTALL	      = cp

LIBRARY	      = libscalar.a

MAKEFILE      = Makefile

OBJS	      = advect.o \
		laplace.o \
		laplace4motion.o \
		scalar.o \
		source.o \
		tdummy.o \
		visco.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = advect.f \
		laplace.f \
		laplace4motion.f \
		scalar.f \
		source.f \
		tdummy.f \
		visco.f

SYSHDRS	      =

#

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
###
advect.o: ../../include/constants.h ../../include/scalar.com \
	../../include/visco.com
scalar.o: ../../include/paramt.h ../../include/constants.h ../../include/three.com \
	../../include/flags.com ../../include/visco.com ../../include/io.com
source.o: ../../include/paramt.h ../../include/constants.h
tdummy.o: ../../include/paramt.h ../../include/three.com ../../include/visco.com
