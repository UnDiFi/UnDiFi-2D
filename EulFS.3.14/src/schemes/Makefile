FFLAGS	      = -I../../include -I.

include $(PETSC_DIR)/lib/petsc/conf/variables
include $(PETSC_DIR)/lib/petsc/conf/rules

DEST	      = $(FSPL_DIR)/lib/$(PETSC_ARCH)/.

EXTHDRS	      = ../../include/flags.com \
		../../include/paramt.h

HDRS	      =

INSTALL	      = cp

LIBRARY	      = libspl.a

MAKEFILE      = Makefile

OBJS	      = petsc_kernel.o \
		limiters.o \
		unsteady1.o \
		unsteady2.o \
		unsteady4.o \
		CDS.o \
		FCT_scheme.o \
		FV_scheme.o \
		FVSys_scheme.o \
		LDA2Sys_scheme.o \
		LDAS_fast.o \
		LDASys_scheme.o \
		LDA_scheme.o \
		LDA2_scheme.o \
		LDASqr_scheme.o \
		LWSys_scheme.o \
		LW_scheme.o \
		LW2_scheme.o \
		LWaSys_scheme.o \
		LWbSys_scheme.o \
		NL2_scheme.o \
		NL_scheme.o \
		NS_scheme.o \
		NSys_scheme.o \
		N_scheme.o \
		PSI_scheme.o \
		SDK_scheme.o \
		SUPG_scheme.o \
		UCVSys_scheme.o \
		VKIB_scheme.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = petsc_kernel.c \
		limiters.f \
		unsteady1.f \
		unsteady2.F \
		unsteady4.F \
		CDS.f \
		FCT_scheme.F \
		FV_scheme.f \
		FVSys_scheme.f \
		LDA2Sys_scheme.F \
		LDAS_fast.F \
		LDASys_scheme.F \
		LDA2_scheme.f \
		LDASqr_scheme.f \
		LWSys_scheme.F \
		LW_scheme.f \
		LW2_scheme.f \
		LWaSys_scheme.F \
		LWbSys_scheme.F \
		NL2_scheme.f \
		NL_scheme.f \
		NS_scheme.f \
		NSys_scheme.f \
		N_scheme.f \
		PSI_scheme.f \
		SDK_scheme.f \
		SUPG_scheme.f \
		UCVSys_scheme.F \
		VKIB_scheme.F

SYSHDRS	      =
#
# I had to explicitely define the rules for CRAY
#
#

all:		$(LIBRARY)

$(LIBRARY):	$(OBJS)
		@echo "Loading $(LIBRARY) ..."
		@ar cru $(LIBRARY) $(OBJS)
		@echo "done"

clobber:;	@rm -f $(OBJS) $(LIBRARY) core tags $(DEST)/$(LIBRARY)

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
CDS.o: ../../include/paramt.h ../../include/flags.com test1.inc test2.inc test0.inc
LDA2Sys_scheme.o: ../../include/paramt.h ../../include/flags.com test1.inc test2.inc test0.inc
LDAS_fast.o: ../../include/paramt.h ../../include/flags.com test1.inc test2.inc test0.inc
LDASys_scheme.o: ../../include/paramt.h ../../include/flags.com test1.inc test2.inc test0.inc
LDA_scheme.o: ../../include/paramt.h
LWSys_scheme.o: ../../include/paramt.h ../../include/flags.com test1.inc test2.inc test0.inc
LWbSys_scheme.o: ../../include/paramt.h ../../include/flags.com test1.inc test2.inc test0.inc
NSys_scheme.o: ../../include/paramt.h ../../include/flags.com test1.inc test2.inc test0.inc
N_scheme.o: ../../include/paramt.h
