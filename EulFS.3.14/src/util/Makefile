include $(PETSC_DIR)/lib/petsc/conf/variables
include $(PETSC_DIR)/lib/petsc/conf/rules

FFLAGS        = -I../../include -c

DEST	      = $(FSPL_DIR)/lib/$(PETSC_ARCH)/.

EXTHDRS	      =

HDRS	      =

INSTALL	      = cp

LIBRARY	      = libutil.a

MAKEFILE      = Makefile

OBJS	      = binsrc.o \
		cross_prod.o \
		csrfmt.o \
		dinit.o \
		div.o \
		dump.o \
		factorial.o \
		gauss.o \
		getunit.o \
		i4rank.o \
		icopy.o \
		icycl.o \
		igthr.o \
		iinit.o \
		isctr.o \
		isdmin.o \
		isortrx.o \
		jcycl.o \
		lenstr.o \
		ludeco.o \
		lusolv.o \
		matins.o \
		matsub.o \
		matsum.o \
		matvec.o \
		mlusol.o \
		mlusolx.o \
		mmio.o \
		nullmat.o \
		poisson.o \
		print.o \
		qsort2.o \
		qsorti.o \
		quadrature.o \
		rnkidx.o \
		r8rank.o \
		reordr.o \
		rg.o \
		simpson.o \
		store.o \
		tecplot.o \
		test.o \
		tet_mesh_order4_adj_count.o \
		toms380.o \
		trapezium.o \
		triangulation_order3_adj_count.o \
		unitmat.o \
		xgemv.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = binsrc.f \
		cross_prod.f \
		csrfmt.f \
		dinit.f \
		div.f \
		dump.f \
		factorial.f \
		gauss.f \
		getunit.f \
		i4rank.f \
		icopy.f \
		icycl.f \
		igthr.f \
		iinit.f \
		isctr.f \
		isdmin.f \
		isortrx.f \
		jcycl.f \
		lenstr.f \
		ludeco.f \
		lusolv.f \
		matins.f \
		matsub.f \
		matsum.f \
		matvec.f \
		mlusol.f \
		mlusolx.f \
		mmio.f \
		nullmat.f \
		poisson.f90 \
		print.f \
		qsort2.f \
		qsorti.f \
		quadrature.f \
		rnkidx.f \
		r8rank.f \
		reordr.f \
		rg.f \
		simpson.f \
		store.f \
		tecplot.f \
		tet_mesh_order4_adj_count.f90 \
		test.f \
		toms380.f \
		trapezium.f \
		triangulation_order3_adj_count.f \
		unitmat.f \
		xgemv.f

SYSHDRS	      =

#
# I had to explicitely define the rules for CRAY
#
#

all:		$(LIBRARY)

$(LIBRARY):	$(OBJS)
		@echo "Loading $(LIBRARY) ..."
		@ar cru $(LIBRARY) $(OBJS)
		@echo "done"

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

checkout:;	@co $(SRCS)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)

