FC            = gfortran

FFLAGS	      = -O

DEST          = ../../bin


EXTHDRS	      =

HDRS	      =

INSTALL	      = cp

LD	      = $(FC)

LDFLAGS	      =

LIBS	      =

MAKEFILE      = Makefile

OBJS	      = main.o \
		seterr.o \
		periodic.o \
		xa112d.o \
		rtri.o

PRINT	      = pr

PROGRAM       = triangle2dat-NEW-$(HOSTTYPE)

SHELL	      = /bin/sh

SRCS	      = main.f \
		seterr.f \
		periodic.f \
		xa112d.f \
		rtri.f

SYSHDRS	      =

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM) -L ../../lib/linux_intel_x86_64 -lport -lfxdr -lmylib  -ltirpc  #-L/opt/intel/mkl721//lib/32 -lmkl_ia32
		@echo "done"

clean:;		@rm -f $(OBJS) core $(PROGRAM)

clobber:;	@rm -f $(OBJS) $(PROGRAM) core tags

depend:;	@mkmf -f $(MAKEFILE)

echo:;		@echo $(HDRS) $(SRCS)

index:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@-strip $(PROGRAM)
		@if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(PROGRAM); $(INSTALL) $(PROGRAM) $(DEST)); fi

print:;		@$(PRINT) $(HDRS) $(SRCS)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)
