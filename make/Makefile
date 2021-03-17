include make.inc

#SRC
SRC      = kinds.f90 errore.f90 constants.f90 global_version.f90 io.f90 date_and_times.f90 \
           randoms.f90 environments.f90 lasercom.f90 parameters.f90 utility.f90 epwcom.f90 \
           readinput.f90 atom_weight.f90 cell_base.f90 ions_base.f90 symm_base.f90 \
           kpoints.f90 wavefct.f90 pwcom.f90 readscf.f90 phcom.f90 readphout.f90 \
           elph2.f90 elphmat.f90 grid.f90  readepw.f90 \
           hamiltonian.f90  surfacehopping.f90 lvcsh.f90            
#           initialsh.f90 dynamics.f90 lvcsh.f90

#objects
OBJ       = $(SRC:.f90=.o)

# ar command and flags - for most architectures: AR = ar, ARFLAGS = ruv
ARCH = ar
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
