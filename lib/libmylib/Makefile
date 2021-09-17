FC       = gfortran
FFLAGS   = -O
DEST     = ../linux_intel_$(HOSTTYPE)/.
EXTHDRS  =
HDRS     =
INSTALL  = cp
LIBRARY  = libmylib.a
MAKEFILE = Makefile
OBJS     = areacoords.o \
					 binsrc.o \
					 binsrch.o \
					 bisrch.o \
					 csrfmt.o \
           dcopy.o \
					 dinit.o \
					 i4rank.o \
					 icopy.o \
					 icycl.o \
					 iinit.o \
					 iperm.o \
					 irand.o \
					 isdmin.o \
					 isortrx.o \
					 jcycl.o \
					 nod2cel.o \
					 nullmat.o \
					 print.o \
					 qsort2.o \
					 qsortd.o \
					 qsorti.o \
					 qsorts.o \
					 r8rank.o \
					 removedups.o \
					 rnkidx.o \
					 rg.o \
					 sortsp.o \
					 neighb2d.o \
					 neighb3d.o \
					 kfact.o \
					 ksolve.o \
					 permat2.o \
					 lenstr.o \
					 gr_line_cross.o \
					 cross_prod.o \
					 area.o \
					 getdia2.o \
					 pspltm.o \
					 MatGetSizeSeq.o \
					 sortn2c.o \
					 toms380.o \
					 tri2renka.o \
					 unitmat.o

PRINT	= pr
SHELL	= /bin/sh
SRCS  = areacoords.f \
				binsrc.f \
				binsrch.f \
				bisrch.f \
				csrfmt.f \
        dcopy.f \
				dinit.f \
				i4rank.f \
				icopy.f \
				icycl.f \
				iinit.f \
				iperm.f \
				irand.f \
				isdmin.f \
				isortrx.f \
				jcycl.f \
				nod2cel.f \
				nullmat.f \
				print.f \
				qsort2.f \
				qsorti.f \
				qsortd.f \
				qsorts.f \
				r8rank.f \
				removedups.f \
				rnkidx.f \
				rg.f \
				sortsp.f \
				neighb2d.f \
				neighb3d.f \
				kfact.f \
				ksolve.f \
				permat2.f \
				lenstr.f \
				gr_line_cross.f \
				cross_prod.f \
				area.f \
				getdia2.f \
				pspltm.f \
				MatGetSizeSeq.f \
				sortn2c.f \
				toms380.f \
				tri2renka.f \
				unitmat.f

SYSHDRS =

all:		$(LIBRARY)

$(LIBRARY):	$(OBJS)
		@echo "Loading $(LIBRARY) ..."
		@ar cru $(LIBRARY) $(OBJS)
		@echo "done"

clean:;		@rm -f $(OBJS) core

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
