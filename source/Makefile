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

OBJS     = fx_msh_sps.o \
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

SRCS    = fx_msh_sps.f \
	        co_utp.f \
	        co_uqp.f \
	        co_urr.f \
		      chg_bnd_ptr.f \
		      main.f \
		      interp.f \
		      rand.f \
		      readmesh.f \
		      readpmap.f \
		      fixendpts.f \
		      rtri.f \
		      fnd_phps.f \
		      rdstrshp.f \
		      setbndrynodeptr.f \
	        re_sdw_info.f \
	        re_inp_data.f \
	        wrt_sdw_info.f \
	        co_norm.f \
	        co_pnt_dipl.f \
					co_state_dps.f \
					fx_state_dps.f \
					mv_dps.f \
	        fx_dps_loc.f \
          ch_sh_topology.f \
					fltr_dls.f \
					rd_dps.f \
					rd_dps_eq.f \
					co_shock.f \
					solzne.f \
					pr_sh_state.f \
					fx_sh_state.f \
	    		wtri.f  \
					wtri0.f \
					wsh_mean.f \
					calc_vel.f \
					mv_grid.f

SYSHDRS	=

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
