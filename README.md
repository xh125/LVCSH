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

## Butterfly of the code

![Alt Butterfly](https://mailustceducn-my.sharepoint.com/:i:/g/personal/xh125_mail_ustc_edu_cn/EWqEBJJC5dxEvMFxXtDiXEoBjtb25NvKd1HZIeMSMUpvJg)

## Recompile Quantum_Espresso

The the code use the EPW output file as main input file. Now only the version of [***V.6.8*** of Quantum_Espresso](https://github.com/QEF/q-e/releases/tag/qe-6.8) be support and need to change the [EPW](https://epw-code.org/) source code and recompile to print out the gmnvkq and vmef in complex formated.

```bash {.line-numbers}
cp LVCSH/docs/QE_change_code/v6.8/* qe-6.8/EPW/src
cd qe-6.8
make epw
```
