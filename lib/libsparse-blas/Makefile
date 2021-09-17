FC            = gfortran

FFLAGS	      = -O

DEST	      = ../linux_intel_x86_64

EXTHDRS	      =

HDRS	      =

INSTALL	      = cp

LIBRARY	      = libsparse-blas.a

MAKEFILE      = Makefile

OBJS	      = compblas.o \
		cp16blas.o \
		doubblas.o \
		realblas.o 

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = compblas.f \
		cp16blas.f \
		doubblas.f \
		realblas.f 

SYSHDRS	      =

all:		$(LIBRARY)

$(LIBRARY):	$(OBJS)
		@echo "Loading $(LIBRARY) ..."
		@ar cru $(LIBRARY) $(OBJS)
		@echo "done"

clean:;		@rm -f $(OBJS) core libsparse-blas.a

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
