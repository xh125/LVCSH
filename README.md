# Fewest Switches Surface Hopping in Materials

A Fortran code to simulation the non-adiabatic Molecular Dynamics in the solid state Materials.
We using a linear vibroninc coupling model implementation of Tully's Fewest Switches Surface Hopping (FSSH) for model problems including
a propagator and an implementation of Tully's model problems described in Tully, [J.C. _J. Chem. Phys._ (1990) **93** 1061](https://aip.scitation.org/doi/abs/10.1063/1.459170).

## Install

The code are writen in Fortran, and need MKL Fortran 95 library. To compiler the code, do as follwing step.

```bash {.line-numbers}
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

```bash {.line-numbers}
cp LVCSH/docs/QE_change_code/v6.8/* qe-6.8/EPW/src
cd qe-6.8
make epw
```
## Example

1. make a work directory  
2. In work dir do epw calculate, the dir structure as follow  
   
```
                    (work directory)  
        ___________________|____________________
       |             |            |             |     
     relax          scf         phonon         epw 
```
In directory epw to calculate the electron-phonon coupling matrix using the changed EPW code. And the output be named dependend on the kpoint: as epw40.out, epw80.out, epw120.out, epw160.out. Used to test the kpoint and qpoint convergence.  
3. 


## Butterfly of the code

![Alt Butterfly](https://github.com/xh125/MarkdownImage/raw/main/Image/Butterfly-lvcsh.png)
