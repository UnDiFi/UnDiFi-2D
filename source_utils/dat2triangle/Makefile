FC            = gfortran

FFLAGS	      = -I./common8/ -O

DEST	      = ../../bin

EXTHDRS	      =

HDRS	      = stack.h

INSTALL	      = cp

LD	      = $(FC)

LDFLAGS	      =

LIBS	      =

MAKEFILE      = Makefile

OBJS	      = main.o \
		readat.o \
		nodcod.o \
		seterr.o \
		solzne.o \
		wtri.o \
	        movep.o \
	        renumber.o \
		wkspce.o

PRINT	      = pr

PROGRAM       = dat2triangle-NEW-$(HOSTTYPE)

SHELL	      = /bin/sh

SRCS	      = main.f \
		readat.f \
		nodcod.f \
		seterr.f \
		solzne.f \
		wtri.f \
		wkspce.f


SYSHDRS	      =

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
		@echo "Linking $(PROGRAM) ..."
		$(LD) $(LDFLAGS) -o $(PROGRAM) $(OBJS) $(LIBS) -L ../../lib/linux_intel_x86_64 -lmylib -lport -lfxdr  -ltirpc

clean:;		@rm -f $(OBJS) core  $(PROGRAM)

clobber:;	@rm -f $(OBJS) $(PROGRAM) core tags

depend:;	@mkmf -f $(MAKEFILE) ROOT=$(ROOT)

echo:;		@echo $(HDRS) $(SRCS)

index:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@-strip $(PROGRAM)
		@if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(PROGRAM); $(INSTALL) $(PROGRAM) $(DEST)); fi

print:;		@$(PRINT) $(HDRS) $(SRCS)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)

$(DEST)/$(PROGRAM): $(SRCS) $(LIBS) $(HDRS) $(EXTHDRS)
		@$(MAKE) -f $(MAKEFILE) ROOT=$(ROOT) DEST=$(DEST) install
