VERSION		= 3.14.0
PROGRAM		= EulFS_$(HOSTTYPE)
all: 		$(PROGRAM)
FFLAGS	        = -I../../include/ -I.
FCPPFLAGS	= $(FCPPFLAGS) -I../../include/ -I.
CFLAGS	        = -I../../include/
FSPL_DIR1=${PWD}
export FSPL_DIR=${FSPL_DIR1}
include ../../../petsc-3.14.6/lib/petsc/conf/variables
include ../../../petsc-3.14.6/lib/petsc/conf/rules

#
# ad hoc fix for SP3
#
#FCPPFLAGS       = ${PETSC_INCLUDE} ${PCONF} ${PETSCFLAGS} ${PETSC_PARCH} \
		${FPPFLAGS} -I$(FSPL_DIR)/include/

DEST		= ../../../bin
INSTALL		= cp
SOURCEC	        =
SOURCEF	        = \
		getidx.f \
		ApplicationFunction.F \
		ApplicationFunction_t.F \
		bndry_iset.F \
		JacobianBoundaryConditions.F \
		RHSFunction.F \
		RHSFunction_t.F \
		RHSJacobian.F \
		RHSJacobian_t.F \
		blockdata.f \
		bndvflx.F \
		clearmem.F \
		lhsbc5.F \
		lhsbc6.F \
		exgeo.F \
		newgeo.F \
		ghost.F \
		ghost2.F \
		init.F \
		iset.F \
		iset_t.F \
		main.F \
		matsch.F \
		MotionSolver.F \
		myTS.F \
		nodres.F \
		nodres_t.F \
		noname.f \
		printmat2.F \
		printmat.F \
		printmatmm.F \
		qb.F \
		rdat.F \
		readat.F \
		rgrdpts.F \
		rhsbc1.F \
		rhsbc4.F \
		rhsbc5.F \
		rhsbc5c.F \
		sclsch.F \
		setbc4laplace.F \
		setibc.F \
		seterr.F \
		setupRHS.F \
		setupRHS_t.F \
		setupLHS_b.F \
		solzne.F \
		MatAllocaSeq.F \
		test.F \
		tmodel.F \
		turbini.F \
		turbsch.F \
		update2.F \
		update3.F \
		update4.F \
		weakbc.F
SOURCEH	        =
OBJSC	        =
OBJSF	        = \
		getidx.o \
		ApplicationFunction.o \
		ApplicationFunction_t.o \
		bndry_iset.o \
		JacobianBoundaryConditions.o \
		RHSFunction.o \
		RHSFunction_t.o \
		RHSJacobian.o \
		RHSJacobian_t.o \
		blockdata.o \
		bndvflx.o \
		clearmem.o \
		lhsbc5.o \
		lhsbc6.o \
		exgeo.o \
		newgeo.o \
		ghost.o \
		ghost2.o \
		init.o \
		iset.o \
		iset_t.o \
		main.o \
		matsch.o \
		MotionSolver.o \
		myTS.o \
		nodres.o \
		nodres_t.o \
		noname.o \
		printmat2.o \
		printmat.o \
		printmatmm.o \
		qb.o \
		rdat.o \
		readat.o \
		rgrdpts.o \
		rhsbc1.o \
		rhsbc4.o \
		rhsbc5.o \
		rhsbc5c.o \
		sclsch.o \
		setbc4laplace.o \
		setibc.o \
		seterr.o \
		setupRHS.o \
		setupRHS_t.o \
		setupLHS_b.o \
		solzne.o \
		MatAllocaSeq.o \
		test.o \
		tmodel.o \
		turbini.o \
		turbsch.o \
		update2.o \
		update3.o \
		update4.o \
		weakbc.o
LIBBASE	        =
#LIBFLAGS      = -L$(HOME)/lib/$(PETSC_ARCH) -lfxdr -lport -lmynag -lskit
LIBFLAGS      = -L ../../../lib/linux_intel_x86_64 -lfxdr -lport -lsparse-blas -lskit -ltirpc
LIBS		= ../../lib/$(PETSC_ARCH)/libscalar.a \
		../../lib/$(PETSC_ARCH)/libeuler.a \
		../../lib/$(PETSC_ARCH)/libspl.a \
		../../lib/$(PETSC_ARCH)/libns.a \
		../../lib/$(PETSC_ARCH)/libturbo.a \
		../../lib/$(PETSC_ARCH)/libgeo.a \
		../../lib/$(PETSC_ARCH)/libchem.a \
		../../lib/$(PETSC_ARCH)/libutil.a
#
# CLDFILES to be defined only for CRAY
#
#CLDFILES	= dp_lapack.cld dp_blas.cld pat.cld
#CLDFILES	= dp_lapack.cld dp_blas.cld

#U_FC         = f77 $(U_FOPTFLAGS)
#U_FOPTFLAGS  = -O

#L_FC         = f77 $(L_FOPTFLAGS)
#L_FOPTFLAGS  = -O

# restore PETSc stuff
#U_FC         = $(FC) $(FOPTFLAGS)
#L_FC         = $(FC) $(FOPTFLAGS)

#FOPTFLAGS    = -O5

$(PROGRAM): $(OBJSF) $(OBJSC) $(LIBS)
	-$(FLINKER) $(CLDFILES) -o $(PROGRAM) $(OBJSF) $(OBJSC) $(LIBS) \
         $(PETSC_FORTRAN_LIB) $(PETSC_LIB) $(LIBFLAGS)


$(FSPL_DIR)/lib/$(PETSC_ARCH)/libgeo.a:
		cd $(FSPL_DIR)/src/geometry; $(MAKE) install
$(FSPL_DIR)/lib/$(PETSC_ARCH)/libeuler.a:
		cd $(FSPL_DIR)/src/euler; $(MAKE) install
$(FSPL_DIR)/lib/$(PETSC_ARCH)/libspl.a:
		cd $(FSPL_DIR)/src/schemes; $(MAKE) install
$(FSPL_DIR)/lib/$(PETSC_ARCH)/libns.a:
		cd $(FSPL_DIR)/src/navier-stokes; $(MAKE) install
$(FSPL_DIR)/lib/$(PETSC_ARCH)/libscalar.a:
		cd $(FSPL_DIR)/src/scalar; $(MAKE) install
$(FSPL_DIR)/lib/$(PETSC_ARCH)/libturbo.a:
		cd $(FSPL_DIR)/src/turbo; $(MAKE) install
$(FSPL_DIR)/lib/$(PETSC_ARCH)/libchem.a:
		cd $(FSPL_DIR)/src/chemistry; $(MAKE) install
$(FSPL_DIR)/lib/$(PETSC_ARCH)/libutil.a:
		cd $(FSPL_DIR)/src/util; $(MAKE) install


#.SUFFIXES:
#.F.o:
#	$(U_FC) -c $(FFLAGS) $(FCPPFLAGS) $<
#.f.o:
#	$(U_FC) -c $(FFLAGS) $<
#########


checkout:;      @co $(SOURCEF)

install:        $(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(PROGRAM); $(INSTALL) $(PROGRAM) $(DEST)); fi
###
blockdata.o: ../../include/paramt.h ../../include/bnd.h \
	../../include/constants.h ../../include/bnd.com \
	../../include/conv.com ../../include/implicit.h \
	../../include/nboun.com ../../include/three.com
bndvflx.o: ../../include/paramt.h ../../include/constants.h \
	../../include/implicit.h ../../include/bnd.h ../../include/bnd.com \
	../../include/three.com ../../include/nloc.com ../../include/flags.com \
	../../include/stream.com ../../include/io.com
exgeo.o: ../../include/io.com ../../include/constants.h ../../include/nloc.com
iset.o: ../../include/iset.com
lhsbc5.o: ../../include/iset.com
main.o:	../../include/stack.com
matsch.o: ../../include/flags.com
mshcnt.o: ../../include/verbose.com ../../include/io.com
nodres.o: ../../include/paramt.h ../../include/bnd.h ../../include/constants.h \
	../../include/bnd.com ../../include/nloc.com ../../include/flags.com \
	../../include/stream.com ../../include/conv.com \
	../../include/nboun.com ../../include/implicit.h ../../include/io.com
psub.o: ../../include/constants.h ../../include/paramt.h ../../include/nloc.com \
	../../include/flags.com
rdat.o: ../../include/paramt.h ../../include/bnd.h ../../include/implicit.h \
	../../include/visco.com ../../include/constants.h \
	../../include/conv.com ../../include/stream.com \
	../../include/chorin.com ../../include/scalar.com \
	../../include/flags.com ../../include/turb.com \
	../../include/bnd.com ../../include/io.com \
	../../include/verbose.com
readat.o: ../../include/constants.h ../../include/bnd.h ../../include/paramt.h \
	../../include/io.com ../../include/nloc.com ../../include/flags.com \
	../../include/stream.com
rhsbc1.o: ../../include/paramt.h ../../include/constants.h \
	../../include/iset.com
rhsbc4.o: ../../include/paramt.h ../../include/constants.h \
	../../include/iset.com
rhsbc5.o: ../../include/paramt.h ../../include/iset.com \
	../../include/constants.h
sclsch.o: ../../include/flags.com
solzne.o: ../../include/io.com
turbcomp.o: ../../include/paramt.h ../../include/constants.h ../../include/nloc.com \
	../../include/three.com ../../include/flags.com ../../include/turb.com \
	../../include/trip.com ../../include/visco.com \
	../../include/nboun.com ../../include/implicit.h ../../include/io.com
update2.o: ../../include/constants.h ../../include/paramt.h \
	../../include/conv.com ../../include/nloc.com ../../include/verbose.com \
	../../include/implicit.h ../../include/iset.com \
	../../include/flags.com ../../include/io.com
update3.o: ../../include/constants.h ../../include/paramt.h \
	../../include/implicit.h ../../include/conv.com ../../include/nloc.com \
	../../include/verbose.com ../../include/iset.com \
	../../include/flags.com ../../include/io.com
update4.o: ../../include/constants.h ../../include/paramt.h \
	../../include/conv.com ../../include/nboun.com ../../include/nloc.com \
	../../include/verbose.com ../../include/implicit.h ../../include/io.com
weakbc.o: ../../include/paramt.h ../../include/constants.h ../../include/bnd.h \
	../../include/bnd.com ../../include/three.com ../../include/nloc.com \
	../../include/implicit.h rotaterhs.f rotaterhs2.f
