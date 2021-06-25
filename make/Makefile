#!/usr/bin/make -f

SHELL = /bin/sh

# path
#project_home = $(dir $(abspath $(MAKEFILE_LIST)))
project_home = ..
srcdir       = $(project_home)/src
bindir       = $(project_home)/bin

# the path to find the source code. 
VPATH = $(srcdir)

# output exe file
exe	    = LVCSH.x

# IFLAGS = how to locate directories with *.h or *.f90 file to be included
IFLAGS         = -I$(project_home)/include

# MOD_FLAG = flag used by f90 compiler to locate modules  
MOD_FLAG      = -I

# Compiler 
F90     = ifort 
MPIF90 	= mpiifort 
CC      = icc 
MPICC   = mpiicc
CXX     = icpc 
MPICXX  = mpiicpc

# C preprocessor and preprocessing flags - for explicit preprocessing, 
# if needed (see the compilation rules above) 
# preprocessing flags must include DFLAGS and IFLAGS
CPP      = cpp
CPPFLAGS = -P -traditional -Uvector $(IFLAGS)

# compiler flags: C, F90
# C flags must include DFLAGS and IFLAGS
# F90 flags must include MODFLAGS, IFLAGS, and FDFLAGS with appropriate syntax
CFLAGS   =  $(OFLAGS) $(IFLAGS)
F90FLAGS =  $(FFLAGS) -fpp $(IFLAGS)
OFLAGS   =  -O3


# compiler flags with and without optimization for fortran-77
# the latter is NEEDED to properly compile dlamch.f, used by lapack
FFLAGS         = -O2
FFLAGS_NOOPT   = -O0 
FDEBUG = -g 
CDEBUG = -g
LDFLAGS= 


# LIB dir
LIBDIR =  
LIBS   = 


#MKL libraries
MKLLIB           = ${MKLROOT}/lib/intel64
MKLINCLUDE       = -I${MKLROOT}/include -I${MKLROOT}/include/intel64/lp64 
BLAS95_LIB       = ${MKLLIB}/libmkl_blas95_lp64.a
LAPACP95_LIB     = ${MKLLIB}/libmkl_lapack95_lp64.a
BLACS_LIBS       = ${MKLLIB}/libmkl_blacs_intelmpi_lp64.a
SCALAPACK_LIBS   = ${MKLLIB}/libmkl_scalapack_lp64.a 
lp64_lib         = ${MKLLIB}/libmkl_intel_lp64.a
thread_lib       = ${MKLLIB}/libmkl_intel_thread.a
sequential_lib   = ${MKLLIB}/libmkl_sequential.a
core_lib         = ${MKLLIB}/libmkl_core.a
MKL_dynamic_Link = ${BLAS95_LIB} ${LAPACP95_LIB} -L${MKLLIB} \
                   -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
                   -liomp5 -lpthread -lm -ldl
MKL_static_Link  = ${BLAS95_LIB} ${LAPACP95_LIB} -Wl,--start-group \
                   ${lp64_lib} ${sequential_lib} ${core_lib} \
									 -Wl,--end-group \
									 -lpthread -lm -ldl

#SRC
SRCS_F= kinds.f90 constants.f90 date_and_times.f90 version.f90 io.f90 environments.f90 \
        pwcom.f90 lasercom.f90 surfacecom.f90 parameters.f90 \
        utility.f90 readinput.f90 \
        cell_base.f90 ions_base.f90 symm_base.f90 readscf.f90 \
        phcom.f90 elph_tetra_mod.f90 readphout.f90 \
        elph2.f90 epwcom.f90 grid.f90 \
        recvec.f90 lrcom.f90 wannierEPW.f90 errore.f90 readepw.f90 \
        randoms.f90  getwcvk.f90 initialsh.f90 \
        hamiltonian.f90 fssh.f90 sc_fssh.f90 cc_fssh.f90 surfacehopping.f90 \
        dynamics.f90 saveinf.f90 lvcsh.f90     					 

SRCS_C=

SRCS = $(SRCS_C) $(SRCS_F)


#---------------------------------------------------------------
#Substitution References
#---------------------------------------------------------------
OBJS       = $(SRCS_F:.f90=.o) $(SRCS_C:.c=.o)


# complilation rules

#--------------------------------------------------------------
# set the suffix list 
#--------------------------------------------------------------
.SUFFIXES:
.SUFFIXES: .mod .o .c .f90 .h .fh

# most fortran compilers can directly preprocess c-like directives: use 
# 	$(MPIF90) $(FFLAGS) -c $< 
# if explicit preprocessing by the C preprocessor is needed, use: 
#	$(CPP) $(CPPFLAGS) $< -o $*.F90 
#	$(MPIF90) $(F90FLAGS) -c $*.F90 -o $*.o 
# remember the tabulator in the first column !!!

#---------------------------------------------------------------
# Redefining Pattern Rules
#---------------------------------------------------------------
%.o : %.f90
	$(F90) $(F90FLAGS) ${MKLINCLUDE} -c $<

%.o : %.c
	$(CC) $(CFLAGS) -c $<

%.fh : %.h
	$(CPP) $(CPPFLAGS) $< -o $*.fh


.PHONY: all

all: $(exe)

##
${exe} : ${OBJS}	bindir
#	${F90} -o ${exe} ${OBJS} ${MKL_dynamic_Link}
	${F90} -o ${exe} ${OBJS} ${MKL_static_Link}
	cp -f ${exe} ../bin

bindir :
	cd ../ ; test -d bin || mkdir bin

#make clean
.PHONY : clean
clean:
	- rm -f ./*.mod
	- rm -f ./*.o
	- rm -rf $(bindir)
	- rm -f ./$(exe)