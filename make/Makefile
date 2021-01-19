include make.inc

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
