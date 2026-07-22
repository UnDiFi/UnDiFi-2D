FC            = gfortran

FFLAGS	      = -O -std=legacy -fno-range-check

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
		rtri.o \
		nodcod.o \
		wtri.o \
		rsu2.o

PRINT	      = pr

PROGRAM       = su22triangle-$(HOSTTYPE)

SHELL	      = /bin/sh

SRCS	      = main.f \
		seterr.f \
		rtri.f \
		nodcod.f \
		wtri.f \
		rsu2.f

SYSHDRS	      =

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM) -L ../../lib/linux_intel_x86_64 -lport -lmylib
		@echo "done"

clean:;		@rm -f $(OBJS) core $(PROGRAM)

clobber:;	@rm -f $(OBJS) $(PROGRAM) core tags

install:	$(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@-strip $(PROGRAM)
		@if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(PROGRAM); $(INSTALL) $(PROGRAM) $(DEST)); fi
