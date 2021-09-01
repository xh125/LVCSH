# Fewest Switches Surface Hopping in Materials

A Fortran code to simulation the non-adiabatic Molecular Dynamics in the Crystal Materials.
We using a linear vibroninc coupling model implementation of Tully's Fewest Switches Surface Hopping (FSSH) for model problems including
a propagator and an implementation of Tully's model problems described in Tully, J.C. _J. Chem. Phys._ (1990) **93** 1061.

## Install

The code are writen in Fortran, and need MKL library. To compiler the code, do as follwing step.

```bash
module load compiler/intel/2020.4.302
tar -zxvf lvcsh.tar.gz
cd LVCSH  
cd make  
make  
cd ../make_complex  
make  
```

## Recompile Quantum_Espresso

The the code use the EPW output file as main input file. Now only the version of [***V.6.8*** of Quantum_Espresso](https://github.com/QEF/q-e/releases/tag/qe-6.8) be support and need to change the [EPW](https://epw-code.org/) source code and recompile to print out the gmnvkq and vmef in complex formated.

```bash
cp LVCSH/docs/QE_change_code/v6.8/* qe-6.8/EPW/src
cd qe-6.8
make epw
```
