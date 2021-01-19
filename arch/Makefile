CMODE = release             # Compiler mode
# options:
# release,
# debug


#compiler   # options:gnu,pgi,intel

COMPILER = intel

# path
PROJECT_HOME = ./
SRCDIR = $(PROJECT_HOME)/src
OBJDIR = $(PROJECT_HOME)/obj_$(COMPILER)/$(CMODE)
MODDIR = $(PROJECT_HOME)/mod
BINDIR = $(PROJECT_HOME)/bin

# PGI fortran and C
ifeq($(COMPILER),pgi)
  FC=mpif90
  LINK=mpif90
  CC=pgcc
  
  ifeq($(CMODE),debug)
    FFLAGS0 = -module $(OBJDIR) -g -r8 -Kieee -C
  else
    CMODE =release
    FFLAGS0 = -module $(OBJDIR) -O2 -r8 -Kieee
  endif
  ARCH =ar
  ARCHFLAG = cr
endif

# Intel fortran and C
ifeq($(COMPILER),intel)
  FC=mpif90
  LINK=mpif90
  CC=icc
  
  ifeq($(CMODE),debug)
    FFLAGS0 = -module $(OBJDIR) -g -r8 -check all
  else
    CMODE = release
    FFLAGS0 = -module $(OBJDIR) -O2 -r8
  endif
  ARCH = ar
  ARCHFLAG = cr
endif

# GNU fortran and C 
ifeq($(COMPILER),gnu)
  FC=mpif90
  LINK=mpif90
  CC=gcc
  
  ifeq(&(CMODE),debug)
    FFLAGS0 =-J$(OBJDIR) -fdefault-real-8 -frecord-marker=4 -g -fcheck=all -fbacktrace
  else
    CMODE=release
    FFLAGS0 = -J$(OBJDIR) -O2 -fdefault-real-8 -frecord-marker=4
  endif
  
  ARCH=ar
  ARCHFLAG = cr
endif


FFLAGS = $(FFLAGS0)


FCC	    =	ifort -fpp # -check all -pg -traceback
FFLAGS      = -g -O2



#MKL libraries
MKLLIB      = -L${MKLROOT}/lib/intel64
MKLINCLUDE  = -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include
FCCFLAG	    = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -liomp5\
              -lmkl_blas95_lp64 -lmkl_lapack95_lp64  


#SRC
SRC      = kinds.f90 constants.f90 control.f90 parameters.f90 io.f90 utility.f90 readinput.f90 \
           softinformation.f90 readposcar.f90 readphonon.f90 readwannierfile.f90 datafitting.f90 \
           hamiltonian.f90 randoms.f90 surfacehopping.f90 initialsh.f90 dynamics.f90 main.f90
#objects
OBJ       = $(SRC:.f90=.o)

# output exe file
PROG	    = excitonSH.x



#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(FCC) $(FFLAGS) ${MKLLIB} ${MKLINCLUDE} ${FCCFLAG} -c $<

##
${PROG} : ${OBJ}	
	${FCC} -o ${PROG} ${OBJ} ${MKLINCLUDE} ${MKLLIB} ${FCCFLAG}
	cp -f ${PROG} ../bin
  
#Compiler utility
Nshift.x : setposcar.o
	${FCC} -o setposcar.o
setposcar.o : setposcar.f90
	${FCC} -c setposcar.f90
	cp Nshift.x ../bin

all:${PROG} Nshift.x
	  
#make clean
.PHONY : clean all
clean:
	-rm ${PROG} ${OBJ} *.mod ../bin/*
	
