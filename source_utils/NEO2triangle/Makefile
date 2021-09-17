FC            = gfortran

FFLAGS	      = -O

DEST	      = ../../bin

EXTHDRS	      =

HDRS	      =          

INSTALL	      = cp

LD	      = $(FC)

LDFLAGS	      =

LIBS	      =

MAKEFILE      = Makefile

OBJS	      = main_d2t.o

PRINT	      = pr

PROGRAM       = NEO2triangle

SHELL	      = /bin/sh

SRCS	      = main_d2t.f90


SYSHDRS	      =

all:		$(PROGRAM)

$(PROGRAM):     $(SRC) $(LIBS)
		@echo "Compling  $(PROGRAM) ..."
		$(FC) $(FFLAGS) $(LDFLAGS) -o $(PROGRAM) $(SRCS) $(LIBS)  

clean:;		@rm -f $(OBJS)  $(PROGRAM)

clobber:;	@rm -f $(OBJS) $(PROGRAM) 

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
