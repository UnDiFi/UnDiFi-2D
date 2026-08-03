#FC      = ifort
FC       = gfortran

#FFLAGS  = -O -Wall -ccdefault fortran
FFLAGS   = -O -Wall

DEST     = ../bin/.

EXTHDRS  =

HDRS     = stack.h

INSTALL  = cp

LD       = $(FC)

LDFLAGS  =

LIBS1    = -lport -lmylib  -lfxdr -lsparse-blas -ltirpc
#LIBS1   = -lport -lmylib  -lfxdr -lmynag -lsparse-blas
#LIBS1   = -lport -lmylib  -lfxdr -lmynag -lsparse-blas -ltirpc

MAKEFILE = Makefile

OBJS     = mod_error.o \
	         fx_msh_sps.o \
	         co_utp.o \
	         co_uqp.o \
	         co_urr.o \
	         chg_bnd_ptr.o \
		       main.o \
		       interp.o \
		       rand.o \
		       readmesh.o \
		       readpmap.o \
		       rtri.o \
		       fnd_phps.o \
		       rdstrshp.o \
		       fixendpts.o \
		       setbndrynodeptr.o \
           re_sdw_info.o \
           re_inp_data.o \
           wrt_sdw_info.o \
           co_norm.o \
           co_pnt_dspl.o \
		       co_state_dps.o \
		       fx_state_dps.o \
		       mv_dps.o \
           fx_dps_loc.o \
		       ch_sh_topology.o \
		       fltr_dls.o \
		       rd_dps.o \
		       rd_dps_eq.o \
		       co_shock.o \
		       solzne.o \
		       pr_sh_state.o \
		       fx_sh_state.o \
		       wtri.o \
		       wtri0.o \
		       wsh_mean.o \
		       calc_vel.o \
		       mv_grid.o

PRINT   = pr

PROGRAM = UnDiFi-2D_$(HOSTTYPE)

SHELL   = /bin/sh

SRCS    = mod_error.f90 \
	        fx_msh_sps.f90 \
	        co_utp.f90 \
	        co_uqp.f90 \
	        co_urr.f90 \
		      chg_bnd_ptr.f90 \
		      main.f90 \
		      interp.f90 \
		      rand.f90 \
		      readmesh.f90 \
		      readpmap.f90 \
		      fixendpts.f90 \
		      rtri.f90 \
		      fnd_phps.f90 \
		      rdstrshp.f90 \
		      setbndrynodeptr.f90 \
	        re_sdw_info.f90 \
	        re_inp_data.f90 \
	        wrt_sdw_info.f90 \
	        co_norm.f90 \
	        co_pnt_dspl.f90 \
					co_state_dps.f90 \
					fx_state_dps.f90 \
					mv_dps.f90 \
	        fx_dps_loc.f90 \
          ch_sh_topology.f90 \
					fltr_dls.f90 \
					rd_dps.f90 \
					rd_dps_eq.f90 \
					co_shock.f90 \
					solzne.f90 \
					pr_sh_state.f90 \
					fx_sh_state.f90 \
	    		wtri.f90  \
					wtri0.f90 \
					wsh_mean.f90 \
					calc_vel.f90 \
					mv_grid.f90

SYSHDRS	=

.SUFFIXES: .f90 .o

.f90.o:
		$(FC) $(FFLAGS) -c $<

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
		@echo "Linking $(PROGRAM) ..."
		$(LD) $(LDFLAGS) -o $(PROGRAM) $(OBJS) $(LIBS1) -L../lib/linux_intel_$(HOSTTYPE)
		@echo "done"

libport.a:
		cd port; $(MAKE) install
libmylib.a:
		cd mylib; $(MAKE) install

clean:;		@rm -f $(OBJS) core *.node *poly *.ele *.neigh *.edge  $(PROGRAM)

clobber:;	@rm -f $(OBJS) $(PROGRAM) core tags *.node *poly *.ele *.neigh *.edge *.dat *.plt $(LIBS); cd port;\
       		make clean; cd ..; cd mylib; make clean

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
