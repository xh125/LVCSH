# CMODE

# options:
# gnu,
# pgi
# intel

COMPILER = gnu

# path
PROJECT_HOME = ..
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

# LIB dir
LIBDIR = #-L$(BLASDIR) -L$(LAPACKDIR)
LIBS = #

# output exe file
EXE=$(BINDIR)/run_$(COMPILER)_$(CMODE)

    