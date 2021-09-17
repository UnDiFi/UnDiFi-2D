FC            = gfortran

FFLAGS	      = -O

DEST	      = ../linux_intel_x86_64

EXTHDRS	      =

HDRS	      =

INSTALL	      = cp

LIBRARY	      = libport.a

MAKEFILE      = Makefile

OBJS	      = d1mach.o \
		e9rint.o \
		entsrc.o \
		eprint.o \
		erroff.o \
		fdump.o \
		i0tk00.o \
		i1mach.o \
		i8save.o \
		istkgt.o \
		istkin.o \
		istkmd.o \
		istkqu.o \
		istkrl.o \
		istkst.o \
		nerror.o \
		r1mach.o \
		retsrc.o \
		s88fmt.o \
		seterr.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = d1mach.f \
		e9rint.f \
		entsrc.f \
		eprint.f \
		erroff.f \
		fdump.f \
		i0tk00.f \
		i1mach.f \
		i8save.f \
		istkgt.f \
		istkin.f \
		istkmd.f \
		istkqu.f \
		istkrl.f \
		istkst.f \
		nerror.f \
		r1mach.f \
		retsrc.f \
		s88fmt.f \
		seterr.f

SYSHDRS	      =

all:		$(LIBRARY)

$(LIBRARY):	$(OBJS)
		@echo "Loading $(LIBRARY) ..."
		@ar cru $(LIBRARY) $(OBJS)
		@echo "done"

clean:;		@rm -f $(OBJS) core libport.a

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
