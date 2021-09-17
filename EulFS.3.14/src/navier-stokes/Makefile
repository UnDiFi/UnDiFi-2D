include $(PETSC_DIR)/lib/petsc/conf/variables
include $(PETSC_DIR)/lib/petsc/conf/rules

FFLAGS	      = -I../../include

DEST	      = $(FSPL_DIR)/lib/$(PETSC_ARCH)/.

EXTHDRS	      = ../../include/bnd.h \
		../../include/constants.h \
		../../include/paramt.h \
		../../include/stream.com \
		../../include/suther.com \
		../../include/three.com \
		../../include/transf.com \
		../../include/turb.com \
		../../include/visco.com

HDRS	      =

INSTALL	      = cp

LIBRARY	      = libns.a

MAKEFILE      = Makefile

OBJS	      = bvsflx2.o \
		bvsflx4.o \
		laminar.o \
		sutherlaw.o \
		vsflx2.o \
		vsflx4.o \
		wskin.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = bvsflx2.f \
		bvsflx4.f \
		laminar.f \
		sutherlaw.f \
		vsflx2.f \
		vsflx4.f \
		wskin.F

SYSHDRS	      =

#
# I had to explicitely define the rules for CRAY
#
#
#.f.o:
#$(FC) $(FFLAGS) $<
#.F.o:
#$(FC) $(FFLAGS) -DPARCH_$(PETSC_ARCH) $<

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
bvsflx2.o: ../../include/paramt.h ../../include/bnd.h ../../include/constants.h \
	../../include/visco.com ../../include/stream.com ../../include/three.com
bvsflx4.o: ../../include/paramt.h ../../include/bnd.h ../../include/constants.h \
	../../include/turb.com ../../include/stream.com \
	../../include/visco.com ../../include/transf.com ../../include/three.com
laminar.o: ../../include/paramt.h ../../include/three.com
sutherlaw.o: ../../include/suther.com
vsflx2.o: ../../include/paramt.h ../../include/constants.h \
	../../include/visco.com ../../include/three.com
vsflx4.o: ../../include/paramt.h ../../include/bnd.h ../../include/constants.h \
	../../include/visco.com ../../include/transf.com \
	../../include/turb.com ../../include/three.com
wskin.o: ../../include/paramt.h ../../include/bnd.h ../../include/bnd.com \
	../../include/stream.com
