# Fewest Switches Surface Hopping

A Fortran code to simulation the non-adiabatic Molecular Dynamics in the Crystal Materials.
We using a linear vibroninc coupling model implementation of Tully's Fewest Switches Surface Hopping (FSSH) for model problems including
a propagator and an implementation of Tully's model problems described in Tully, J.C. _J. Chem. Phys._ (1990) **93** 1061.

## Install

1. tar -zxvf lvcsh.tar.gz  
2. module load compiler/intel/2020.4.302  
3. cd LVCSH  
4. cd make  
5. make  
6. cd ../make_complex  
7. make  

## Recomple Quantum_Espresso

The the code use the EPW output file as main input file. Now only the version of ***V.6.8*** of Quantum_Espresso be support and need to change the EPW source code and recompile to print out the gmnvkq and vmef in complex formated. 

1. cp LVCSH/docs/QE_change_code/v6.8/* qe-6.8/EPW/src
2. cd qe-6.8
3. make epw

