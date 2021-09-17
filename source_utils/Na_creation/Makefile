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

OBJS	      = physical_values.mod \
		reading_mod.mod \
		writing_mod.mod

PRINT	      = pr

PROGRAM       = Na_creation

SHELL	      = /bin/sh

SRCS	      = physical_values.f90 \
		reading_mod.f90 \
		writing_mod.f90 \
		main.f90 

SYSHDRS	      =

all:		$(PROGRAM)

$(PROGRAM):     $(SRCS)  $(LIBS)
		@echo "Compiling $(PROGRAM) ..."
		$(FC) $(FFLAGS) -o $(PROGRAM) $(SRCS) $(LIBS)  

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
