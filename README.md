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

>In directory epw to calculate the electron-phonon coupling matrix using the changed EPW code. And the output be named dependend on the kpoint: as epw40.out, epw80.out, epw120.out, epw160.out. Used to test the kpoint and qpoint convergence.  

3. make a directory for lvcsh calculation  

```bash
mkdir LVCSH
```

4. make a lsf job script. Need to change the BUSB -q,-n,-R and MODULEPATH as your Environment.   

```bash
lvcsh.bsub
#!/bin/bash
#BSUB -J lvcsh-epw
#BSUB -q privateq-zw
#BSUB -n 32
#BSUB -R "span[ptile=32]"
#BSUB -o %J.out
#BSUB -e %J.err

export MODULEPATH=/share/home/zw/xiehua/opt/modules-4.7.1/modulefiles

module load lvcsh/0.6.0

CURDIR=$PWD
#Generate nodelist
rm -f ${CURDIR}/nodelist >& /dev/null
for i in `echo $LSB_HOSTS`
do
echo $i >> ${CURDIR}/nodelist
done

NP=`cat ${CURDIR}/nodelist |wc -l`

for i in {1..${NP}}
do
mkdir sample$i
cp ./LVCSH.in ./sample$i/
cd ./sample$i
LVCSH_complex.x &
cd ..
done
wait

```

5. Use shell script mkepwdir.sh to build dir for diffrent kpoints and qpoints Surface hopping calculation. And the script will make a test running in the different director to give how to set nefre_sh and nhfre_sh.  

```bash
mkepwdir.sh
#!/bin/bash
ncore=28
for i in $(seq 40 40 160)
    do 
        mkdir epw$i
        mkdir epw$i/QEfiles
        cp ../epw/epw$i.out epw$i/QEfiles/
        cp lvcsh.bsub epw$i
        sed -i "s/ncore/$ncore/g" epw$i/lvcsh.bsub
        sed -i "2s/lvcsh-epw/lvcsh-epw$i-n0/g" epw$i/lvcsh.bsub
        cp job.sh epw$i
        cp LVCSH.in epw$i
        sed -i "s:epw40:epw$i:g" epw$i/LVCSH.in
        cp LVCSH.in epw$i/QEfiles
        sed -i "s/low/high/g" epw$i/QEfiles/LVCSH.in
        sed -i "s:../../QEfiles/epw40.out:epw$i.out:g" epw$i/QEfiles/LVCSH.in
        cp lvcsh-test.bsub epw$i/QEfiles
        cd epw$i/QEfiles
        sed -i "2s/lvcsh-epw40/lvcsh-epw$i-test/g" lvcsh-test.bsub
        bsub < lvcsh-test.bsub
        cd ../..
    done
```  

```bash
job.sh
#!/bin/bash
    for i in {1..10}
    do 
    mkdir node$i
    cp ./lvcsh.bsub ./node$i/
    sed -i "2s/n0/n$i/g" ./node$i/lvcsh.bsub
    cp ./LVCSH.in ./node$i/
    cd ./node$i
    bsub < lvcsh.bsub
    cd ..
    done  
```

```fortran
LVCSH.in
calculation   = "lvcsh"
verbosity     = "low"
outdir        = "./"
methodsh      = "FSSH"
ldecoherence  = .true.
Cdecoherence  = 0.1
lit_gmnvkq    = 0.0    ! in unit of meV
lit_ephonon   = 1.0    ! in unit of meV
lfeedback     = .true.
lehpairsh     = .true.
!lelecsh       = .true.
!lholesh       = .true.
!ieband_min    = 3
!ieband_max    = 4
!ihband_min    = 1
!ihband_max    = 2
!lsortpes      = .false.
!mix_thr       = 0.8
!nefre_sh      = 40
!nhfre_sh      = 40
epwoutname    = "./QEfiles/epw40.out"
naver         = 10
nstep         = 2
nsnap         = 2
dt            = 0.5
pre_nstep     = 50000
pre_dt        = 0.5
gamma         = 0.0   ! in unit of ps-1
ld_fric       = 0.01 !
temp          = 300
l_ph_quantum  = .true.
llaser        = .true.
efield_cart   = 1.0 1.0 1.0
w_laser       = 2.0  ! in unit of eV
fwhm          = 100  ! in unit of fs
nnode         = 1
ncore         = 32
savedsnap     = 1
```  

```bash
lvcsh-test.bsub
#!/bin/bash
#BSUB -J lvcsh-epw40
#BSUB -q privateq-zw
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -o %J.out
#BSUB -e %J.err

export MODULEPATH=/share/home/zw/xiehua/opt/modules-4.7.1/modulefiles
module load lvcsh/0.6.2

LVCSH_complex.x

```  

6. By look the initial adiabatic state in the QEfiles/LVCSH.out for different kpoints directory. Set the **nefre_sh** and **nhfre_sh** in the QEfiles/LVCSH.in to tests the time for one step nonadiabatic calculation. Then, subscrib the job again.  

```shell
bsub < lvcsh-test.bsub
```

7. 

## Butterfly of the code

![Alt Butterfly](https://github.com/xh125/MarkdownImage/raw/main/Image/Butterfly-lvcsh.png)
