FFLAGS	      = -I../../include

include $(PETSC_DIR)/lib/petsc/conf/variables
include $(PETSC_DIR)/lib/petsc/conf/rules

DEST	      = $(FSPL_DIR)/lib/$(PETSC_ARCH)

EXTHDRS	      = ../../include/blkEulerII.com \
		../../include/bnd.h \
		../../include/chorin.com \
		../../include/constants.h \
		../../include/conv.com \
		../../include/flags.com \
		../../include/frame.com \
		../../include/implicit.h \
		../../include/nloc.com \
		../../include/paramt.h \
		../../include/stream.com \
		../../include/three.com \
		../../include/transf.com

HDRS	      =

INSTALL	      = cp

LIBRARY	      = libeuler.a

MAKEFILE      = Makefile

OBJS	      = EulerII.o \
		EulerIIbis.o \
		EulerVII.o \
		EulerVIII.o \
		EulerIX.o \
		EulerX.o \
		EulerXI.o \
		EulerXII.o \
		MatrixSplitting.o \
		bc28vii.o \
		bc28vii4Ar.o \
		bc28viii.o \
		bc4vii.o \
		bc5vii.o \
		bc5vii4Ar.o \
		bc5viii.o \
		beta.o \
		check.o \
		cons2parm.o \
		cons2parm4Ar.o \
		cons2symm.o \
		eigen_II.o \
		eigen_IIbis.o \
		eigen_X.o \
		eigen_IX.o \
		eigen_VII.o \
		eigen_VIII.o \
		eigen_XI.o \
		eigen_XII.o \
		flxale.o \
		flxb4.o \
		flxw2.o \
		flxw4.o \
		frame.o \
		gradprim.o \
		linavg.o \
		linearize.o \
		matsplII.o \
		matsplVII.o \
		matsplVIII.o \
		matsplIX.o \
		matsplX.o \
		matsplXI.o \
		matsplXII.o \
		parm2charHE.o \
		parm2cons.o \
		parm2cons4Ar.o \
		parm2prim.o \
		parm2symm.o \
		pforce.o \
		pressc.o \
		pressi.o \
		pressderiv.o \
		pressderiv.o \
		symm2cons.o \
		symm2parm.o \
		transf.o \
		vstate.o \
		wall.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = EulerII.F \
		EulerIIbis.f \
		EulerIX.f \
		EulerVII.f \
		EulerVIII.f \
		EulerX.f \
		EulerXI.f \
		EulerXII.f \
		MatrixSplitting.f \
		bc28vii.F \
		bc28vii4Ar.F \
		bc28viii.F \
		bc4vii.F \
		bc5vii.F \
		bc5vii4Ar.F \
		bc5viii.f \
		beta.f \
		check.f \
		cons2parm.f \
		cons2parm4Ar.f \
		cons2symm.f \
		eigen_II.F \
		eigen_IIbis.F \
		eigen_IX.F \
		eigen_VII.F \
		eigen_VIII.f \
		eigen_X.F \
		eigen_XI.F \
		eigen_XII.F \
		flxale.f \
		flxb4.f \
		flxw2.f \
		flxw4.f \
		frame.f \
		gradprim.f \
		linavg.f \
		linearize.f \
		matsplII.F \
		matsplVII.F \
		matsplVIII.F \
		matsplIX.F \
		matsplX.F \
		matsplXI.F \
		matsplXII.F \
		parm2charHE.f \
		parm2cons.f \
		parm2cons4Ar.f \
		parm2prim.f \
		parm2symm.f \
		pforce.f \
		pressc.f \
		pressi.f \
		pressderiv.f \
		symm2cons.f \
		symm2parm.f \
		transf.f \
		vstate.f \
		wall.f

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
EulerII.o: ../../include/paramt.h ../../include/bnd.h ../../include/constants.h \
	../../include/blkEulerII.com ../../include/three.com ../../include/transf.com \
	../../include/flags.com
EulerIIbis.o: ../../include/paramt.h ../../include/constants.h \
	../../include/bnd.h ../../include/blkEulerII.com ../../include/three.com \
	../../include/transf.com ../../include/flags.com
EulerIX.o: ../../include/paramt.h ../../include/constants.h ../../include/bnd.h \
	../../include/three.com ../../include/transf.com ../../include/flags.com
EulerVII.o: ../../include/paramt.h ../../include/constants.h \
	../../include/bnd.h ../../include/three.com ../../include/transf.com \
	../../include/flags.com
EulerVIII.o: ../../include/paramt.h ../../include/constants.h \
	../../include/bnd.h ../../include/three.com ../../include/transf.com \
	../../include/flags.com
MatrixSplitting.o: ../../include/paramt.h ../../include/flags.com
bc5vii.o: ../../include/paramt.h ../../include/constants.h \
	../../include/implicit.h ../../include/bnd.h ../../include/stream.com \
	../../include/three.com
bc5viii.o: ../../include/paramt.h ../../include/constants.h ../../include/bnd.h \
	../../include/stream.com
check.o: ../../include/paramt.h ../../include/constants.h ../../include/three.com \
	../../include/chorin.com ../../include/flags.com
cons2parm.o: ../../include/constants.h
eigen_II.o: ../../include/paramt.h ../../include/constants.h \
	../../include/three.com ../../include/blkEulerII.com ../../include/frame.com
eigen_IIbis.o: ../../include/paramt.h ../../include/constants.h \
	../../include/blkEulerII.com ../../include/frame.com ../../include/three.com
eigen_IX.o: ../../include/paramt.h ../../include/constants.h \
	../../include/chorin.com ../../include/three.com ../../include/frame.com
eigen_VII.o: ../../include/paramt.h ../../include/constants.h \
	../../include/three.com
eigen_VIII.o: ../../include/paramt.h ../../include/constants.h \
	../../include/three.com ../../include/chorin.com
frame.o: ../../include/paramt.h ../../include/constants.h ../../include/three.com \
	../../include/frame.com
linearize.o: ../../include/paramt.h ../../include/constants.h \
	../../include/three.com
matsplII.o: ../../include/paramt.h ../../include/blkEulerII.com \
	../../include/frame.com ../../include/constants.h
matsplIX.o: ../../include/paramt.h ../../include/constants.h \
	../../include/frame.com
matsplVII.o: ../../include/paramt.h ../../include/three.com \
	../../include/flags.com ../../include/constants.h
matsplVIII.o: ../../include/paramt.h ../../include/three.com \
	../../include/chorin.com ../../include/flags.com \
	../../include/constants.h
pforce.o: ../../include/paramt.h ../../include/bnd.h ../../include/stream.com
pressc.o: ../../include/constants.h
transf.o: ../../include/constants.h ../../include/paramt.h ../../include/three.com
