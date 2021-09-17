include $(PETSC_DIR)/lib/petsc/conf/variables
include $(PETSC_DIR)/lib/petsc/conf/rules

FFLAGS	      = -I../../include -c

DEST	      = $(FSPL_DIR)/lib/$(PETSC_ARCH)/.

EXTHDRS	      = ../../include/bnd.com \
		../../include/bnd.h \
		../../include/constants.h \
		../../include/io.com \
		../../include/nloc.com

HDRS	      =

INSTALL	      = cp

LIBRARY	      = libgeo.a

MAKEFILE      = Makefile

OBJS	      = aratio.o \
		avgflx.o \
		backup.o \
		bndchk.o \
		cdist.o \
		chkflx.o \
		chkflx2.o \
		chkgcl.o \
		chkint.o \
		celptr.o \
		ff.o \
		getvol.o \
		MatGetSizeSeq.o \
		median.o \
		motion_type.o \
		mshcnt.o \
		nodcod.o \
		periodic.o \
		probe.o \
		rankcl.o \
		rotate.o \
		seg3d.o \
		setbndflxptr.o \
		subyy.o \
		subyz.o \
		setpto.o \
		updtp.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = aratio.f \
		avgflx.f \
		backup.f \
		bndchk.f \
		cdist.f \
		chkflx.f \
		chkflx2.f \
		chkgcl.f \
		chkint.f \
		celptr.f \
		ff.f \
		getvol.f \
		MatGetSizeSeq.f \
		median.f \
		motion_type.f \
		mshcnt.F \
		nodcod.f \
		periodic.F \
		probe.f \
		rankcl.f \
		rotate.f \
		seg3d.f \
		setbndflxptr.f \
		subyy.f \
		subyz.f \
		setpto.f \
		updtp.f

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

checkout:;	@co $(SRCS)

print:;		@$(PRINT) $(HDRS) $(SRCS)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)
###
backup.o: ../../include/nloc.com
bndchk.o: ../../include/bnd.h ../../include/constants.h ../../include/bnd.com \
	../../include/io.com
ff.o: ../../include/constants.h ../../include/io.com
getvol.o: ../../include/constants.h ../../include/io.com
nodcod.o: ../../include/io.com
setpto.o: ../../include/bnd.h ../../include/bnd.com ../../include/io.com
