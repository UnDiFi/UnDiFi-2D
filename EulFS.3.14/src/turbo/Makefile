FFLAGS	      = -I../../include

include $(PETSC_DIR)/lib/petsc/conf/variables
include $(PETSC_DIR)/lib/petsc/conf/rules

DEST	      = $(FSPL_DIR)/lib/$(PETSC_ARCH)/.

EXTHDRS	      = ../../include/bnd.com \
		../../include/bnd.h \
		../../include/constants.h \
		../../include/flags.com \
		../../include/nboun.com \
		../../include/nloc.com \
		../../include/paramt.h \
		../../include/stream.com \
		../../include/three.com \
		../../include/turb.com \
		../../include/visco.com

HDRS	      =

INSTALL	      = cp

LIBRARY	      = libturbo.a

MAKEFILE      = Makefile

OBJS	      = limitvisct.o \
		sa.o \
		sa7.o \
		sa8.o \
		sa_neg_Bassi.o \
		tbdest.o \
		tbdiff.o \
		tbprod.o \
		tbtrip.o \
		tdiff.o \
		tdu.o \
		tft1.o \
		tft2.o \
		tfv1.o \
		tfv2.o \
		tfv3.o \
		tfw.o \
		tgt.o \
		tsa.o \
		tsa1.o \
		tsc.o \
		tsetup.o \
		turbcomp.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = limitvisct.f \
		sa.F \
		sa7.F \
		sa8.F \
		sa_neg_Bassi.F \
		tbdest.f \
		tbdiff.f \
		tbprod.f \
		tbtrip.f \
		tdiff.f \
		tdu.f \
		tft1.f \
		tft2.f \
		tfv1.f \
		tfv2.f \
		tfv3.f \
		tfw.f \
		tgt.f \
		tsa.f \
		tsa1.f \
		tsc.f \
		tsetup.f \
		turbcomp.F

SYSHDRS	      =

#
# I had to explicitely define the rules for CRAY
#
#
#.f.o:
#		$(FC) $(FFLAGS) $<
#.F.o:
#		$(FC) $(FFLAGS) -DPARCH_$(PETSC_ARCH) $<

all:		$(LIBRARY)

$(LIBRARY):	$(OBJS)
		@echo "Loading $(LIBRARY) ..."
		@ar cru $(LIBRARY) $(OBJS)
		@echo "done"

clobber:;	@rm -f $(OBJS) $(LIBRARY) core tags log

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
tbdest.o: ../../include/turb.com ../../include/visco.com
tbdiff.o: ../../include/paramt.h ../../include/three.com ../../include/turb.com \
	../../include/visco.com
tbprod.o: ../../include/turb.com
tbtrip.o: ../../include/visco.com
tdiff.o: ../../include/paramt.h ../../include/three.com
tft1.o: ../../include/turb.com
tft2.o: ../../include/turb.com
tfv1.o: ../../include/turb.com
tfv2.o: ../../include/turb.com ../../include/constants.h
tfw.o: ../../include/turb.com
tgt.o: ../../include/turb.com
tindx.o: ../../include/paramt.h ../../include/bnd.h ../../include/bnd.com \
	../../include/nloc.com ../../include/three.com ../../include/stream.com \
	../../include/turb.com ../../include/visco.com
tsa1.o: ../../include/paramt.h ../../include/three.com
tsc.o: ../../include/paramt.h ../../include/bnd.h ../../include/turb.com \
	../../include/stream.com ../../include/visco.com
tsetup.o: ../../include/paramt.h ../../include/bnd.h ../../include/three.com \
	../../include/stream.com ../../include/turb.com ../../include/nloc.com
turbcomp.o: ../../include/paramt.h ../../include/constants.h ../../include/nloc.com \
	../../include/three.com ../../include/flags.com ../../include/turb.com \
	../../include/visco.com ../../include/nboun.com
