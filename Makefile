# Compiler mode
CMODE = release             
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

# output exe file
EXE	    = $(BINDIR)/LVCSH_$(COMPILER)_$(CMODE).x

# PGI fortran and C
ifeq($(COMPILER),pgi)
  FC  = mpif90
  CC  = pgcc
  CXX = 
  
  ifeq($(CMODE),debug)
    FFLAGS0 = -module $(OBJDIR) -g -r8 -Kieee -C
  else
    CMODE =release
    FFLAGS0 = -module $(OBJDIR) -O2 -r8 -Kieee
  endif

endif

# Intel fortran and C
ifeq($(COMPILER),intel)
  FC   =mpiifort
  CC   =icc
  CXX  =icpc
  
  ifeq($(CMODE),debug)
    FFLAGS0 = -module $(OBJDIR) -g -r8 -check all
  else
    CMODE = release
    FFLAGS0 = -module $(OBJDIR) -O2 -r8
  endif

endif

# GNU fortran and C 
ifeq($(COMPILER),gnu)
  FC  = mpif90
  CC  = gcc
  CXX =
  
  ifeq(&(CMODE),debug)
    FFLAGS0 =-J$(OBJDIR) -fdefault-real-8 -frecord-marker=4 -g -fcheck=all -fbacktrace
  else
    CMODE=release
    FFLAGS0 = -J$(OBJDIR) -O2 -fdefault-real-8 -frecord-marker=4
  endif
  
endif


FFLAGS = $(FFLAGS0)
CFLAGS = 
OFLAG  = -O2
DEBUG  = -O0

# LIB dir 这里加入必要的库和路径
LIBDIR= #-L$(BLASDIR) -L$(LAPACKDIR) 
LIBS   =  # 

#MKL libraries
MKL_PATH    = ${MKLROOT}/lib/intel64
BLAS        =
LAPACK      =
BLACS       = -lmkl_blacs_intelmpi_lp64
SCALAPACK   = $(MKL_PATH)/libmkl_scalapack_lp64.a $(BLACS)
MKL_INCLUDE = -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include
FCFLAG	    = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -liomp5\
              -lmkl_blas95_lp64 -lmkl_lapack95_lp64  

GENCODE_ARCH    := -gencode=arch=compute_30,code=\"sm_30,compute_30\" \
                   -gencode=arch=compute_35,code=\"sm_35,compute_35\" \
                   -gencode=arch=compute_60,code=\"sm_60,compute_60\" \
                   -gencode=arch=compute_70,code=\"sm_70,compute_70\" \
                   -gencode=arch=compute_72,code=\"sm_72,compute_72\"

MPI_INC    = $(I_MPI_ROOT)/include64/


#SRC
SRC      = kinds.f90 constants.f90 control.f90 parameters.f90 io.f90 utility.f90 readinput.f90 \
           softinformation.f90 readposcar.f90 readphonon.f90 readwannierfile.f90 datafitting.f90 \
           hamiltonian.f90 randoms.f90 surfacehopping.f90 initialsh.f90 dynamics.f90 main.f90

#objects
OBJ       = $(SRC:.f90=.o)

# ar command and flags - for most architectures: AR = ar, ARFLAGS = ruv
ARCH =ar
ARCHFLAG = ruv

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(FCC) $(FFLAGS) ${MKLLIB} ${MKLINCLUDE} ${FCCFLAG} -c $<

##
${EXE} : ${OBJ}	
	${FC} -o ${EXE} ${OBJ} ${MKLINCLUDE} ${MKLLIB} ${FCCFLAG}

all:${EXE}
	  
#make clean
.PHONY : clean all
clean:
	rm -f $(OBJDIR)/*.mod
	rm -f $(OBJDIR)/*.o
	rm -f $(BINDIR)/*

#make object dir
dir:
	mkdir ./obj_$(COMPILER)
	mkdir ./obj_$(COMPILER)/debug
	mkdir ./obj_$(COMPILER)/release
