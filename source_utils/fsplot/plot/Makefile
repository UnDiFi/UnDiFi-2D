FC            = gfortran

FFLAGS	      = -I./common8/ -Wall -mcmodel=large
FFLAGS	      = -I./common8/ -Wall

DEST	      = ../../../bin/

EXTHDRS	      = common8/IO \
		common8/bnd \
		common8/bnd.h \
		common8/constants \
		common8/dim_flags \
		common8/es.h \
		common8/int_flags \
		common8/mesh_i4

HDRS	      = stack.h

INSTALL	      = cp

LD	      = $(FC)

LDFLAGS	      =

LIBS	      =  $(TEC_360_2008)/lib/tecio64.a $(TEC_360_2008)/lib/mesa/libGL.so.1 /usr/lib64/libXm.so.4
LIBS	      =  $(TEC_360_2008)/lib/tecio64.a $(TEC_360_2008)/lib/mesa/libtec.so $(TEC_360_2008)/lib/mesa/libGL.so.1 $(TEC_360_2008)/lib/mesa/libOSMesa.so.5 /usr/lib64/libXm.so.4 $(TEC_360_2008)/lib/libgltt.so.0 $(TEC_360_2008)/lib/tcl8.4/libtcl8.4.so $(TEC_360_2008)/lib/tcl8.4/libtk8.4.so $(TEC_360_2008)/lib/libttf.so.2
LIBS	      =  $(TEC_360_2008)/lib/tecio64.a $(TEC_360_2008)/lib/libtec.so $(TEC_360_2008)/lib/libgltt.so.0 /usr/lib64/libXm.so.4 $(TEC_360_2008)/lib/tcl8.4/libtcl8.4.so $(TEC_360_2008)/lib/tcl8.4/libtk8.4.so $(TEC_360_2008)/lib/libttf.so.2
LIBS	      =  $(TEC_360_2008)/lib/tecio64.a $(TEC_360_2008)/lib/libtec.so $(TEC_360_2008)/lib/libgltt.so.0 /usr/lib64/libXm.so.4 $(TEC_360_2008)/lib/tcl8.4/libtcl8.4.so $(TEC_360_2008)/lib/tcl8.4/libtk8.4.so $(TEC_360_2008)/lib/libttf.so.2 /usr/lib64/libGLU.so  $(TEC_360_2008)/lib/mesa/libGL.so.1 $(TEC_360_2008)/lib/mesa/libOSMesa.so.5
LIBS	      = 

MAKEFILE      = Makefile

OBJS	      = ascii2d.o \
		avs.o \
		binsrch2.o \
		bndplot.o \
		bl.o \
		Burkardt.o \
		cf.o \
		ch3d.o \
		cut.o \
		cut3d.o \
		dplot.o \
		gmv.o \
		hydra.o \
		hdb.o \
		interp.o \
		intp3d.o \
		inversed.o \
		main.o \
		mshsiz.o \
		nodcod.o \
		plane.o \
		plpt.o \
		xaa23s.o \
		readat.o \
		rotate2d.o \
		rsub.o \
		setbndrypntr.o \
		setbndrynodepntr.o \
		seterr.o \
		solzne.o \
		sutherlaw.o \
		tecplot.o \
		transf.o \
		wascii.o \
		wkspce.o \
		wmesh.o \
		wtec2.o \
		ydist.o \
		xnufac.o

PRINT	      = pr

PROGRAM       = fsplot-$(HOSTTYPE)

SHELL	      = /bin/sh

SRCS	      = ascii2d.f \
		avs.f \
		binsrch2.f \
		bndplot.f \
		bl.f \
		Burkardt.f \
		cf.f \
		ch3d.f \
		cut.f \
		cut3d.f \
		dplot.f \
		gmv.f \
		hydra.f \
		hdb.f \
		interp.f \
		intp3d.f \
		inversed.f \
		main.f \
		mshsiz.f \
		nodcod.f \
		plane.f \
		plpt.f \
		xaa23s.f \
		readat.f \
		rotate2d.f \
		rsub.f \
		setbndrypntr.f \
		setbndrynodepntr.f \
		seterr.f \
		solzne.f \
		sutherlaw.f \
		tecplot.f \
		transf.f \
		wascii.f \
		wkspce.f \
		wmesh.f \
		wtec2.f \
		ydist.f \
		xnufac.f

SYSHDRS	      =

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
		@echo "Linking $(PROGRAM) ..."
		$(LD) $(LDFLAGS) -o $(PROGRAM) $(OBJS) $(LIBS) -L/usr/local/lib -L ../../../lib/linux_intel_x86_64 -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lm -lfxdr -ltirpc -lport -lfblas -lmylib -lgeometryf90 -ltriangulationf90 -ltriangulation 
		@echo "done"

clean:;		@rm -f $(OBJS) core

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
###
avs.o: common8/constants common8/mesh_i4 common8/int_flags common8/dim_flags \
	common8/IO
blockdata.o: common8/constants common8/bnd.h common8/bnd
bndplot.o: common8/constants common8/bnd.h common8/bnd common8/mesh_i4 \
	common8/int_flags common8/dim_flags common8/IO
dplot.o: common8/constants common8/bnd.h common8/bnd common8/dim_flags \
	common8/int_flags common8/mesh_i4 common8/IO
main.o: common8/constants common8/IO common8/es.h common8/mesh_i4 \
	common8/int_flags common8/dim_flags
readat.o: common8/es.h common8/IO common8/dim_flags common8/int_flags \
	common8/mesh_i4
tecplot.o: common8/constants common8/mesh_i4 common8/int_flags \
	common8/dim_flags common8/IO
transf.o: common8/constants common8/dim_flags common8/int_flags \
	common8/mesh_i4
wkspce.o: common8/mesh_i4 common8/int_flags common8/dim_flags common8/IO \
	common8/es.h
