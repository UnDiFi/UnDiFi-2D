# $Id: makefile,v 1.19 2001/04/10 19:35:58 bsmith Exp $ 



CFLAGS 	        = 
FFLAGS	        = 
CPPFLAGS        =
FPPFLAGS        =
LOCDIR          = src/mat/examples/tutorials/
EXAMPLESC       = ex1.c
EXAMPLESF       =
MANSEC          = Mat

include ${PETSC_DIR}/bmake/common/base

ex1: ex1.o  chkopts
	-${CLINKER}  -o ex1 ex1.o ${PETSC_MAT_LIB} 
	${RM} ex1.o

ex2: ex2.o icopy.o iinit.o chkopts
	-${CLINKER}  -o ex2 ex2.o icopy.o iinit.o ${PETSC_MAT_LIB} /home/aldo/src/petsc-2.3.0/lib/linux_gnu/libgpsk.a
	${RM} ex2.o

#-------------------------------------------------------------
TESTEXAMPLES_1  = ex1.PETSc ex1.rm
TESTEXAMPLES_2  =
TESTEXAMPLES_3  =
TESTEXAMPLES_4  =
TESTEXAMPLES_6  = 
TESTEXAMPLES_7  =
TESTEXAMPLES_9  =
TESTEXAMPLES_13 = ex1.PETSc ex1.rm


include ${PETSC_DIR}/bmake/common/test

