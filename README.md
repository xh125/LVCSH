# Fewest Switches Surface Hopping in Materials

A Fortran code to simulation the non-adiabatic Molecular Dynamics in the solid state Materials.
We using a linear vibroninc coupling model implementation of Tully's Fewest Switches Surface Hopping (FSSH) for model problems including
a propagator and an implementation of Tully's model problems described in Tully, [J.C. _J. Chem. Phys._ (1990) **93** 1061](https://aip.scitation.org/doi/abs/10.1063/1.459170).

## Install

The code are writen in Fortran, and need MKL Fortran 95 library. To compiler the code, do as follwing step.

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

The the code use the EPW output file as main input file. Now only the version of [**V.6.8** of Quantum_Espresso](https://github.com/QEF/q-e/releases/tag/qe-6.8) be support and need to change the [EPW](https://epw-code.org/) source code and recompile to print out the gmnvkq and vmef in complex formated.

```bash
cp LVCSH/docs/QE_change_code/v6.8/* qe-6.8/EPW/src
cd qe-6.8
make epw
```  

## Tutorial

1. make a work directory  

2. In work directory do epw calculate, to get the electron-phonon coupling matrix. the dir structure as follow  

   ```bash
                    (work directory)  
        ___________________|_______________________________
       |          |     |           |             |        |
     relax       scf   dos        phonon         epw      LVCSH
   ```  

   计算电声耦合项的步骤如下所示。在不同的目录下面进行结构优化、自洽计算、声子谱计算和电声耦合计算  

   2.1 进入relax目录，创建vc-relax.in 进行晶格结构优化。对于参数的设置，`etot_conv_thr`,`forc_conv_thr`,`conv_thr`和`press_conv_thr`的精度需要设置的比较高，需要用于后续计算**phonon**. 在该案例中收敛判据可能不够，在实际的计算中，应该根据需要进行设置。
  
   ```fortran
   vc-relax.in
   &CONTROL
    calculation   = "vc-relax"  
    restart_mode  = "from_scratch"
    prefix        = "carbyne"
    outdir        = "./outdir/"
    pseudo_dir    = "./pseudo/"
    verbosity     = "high"
    tprnfor       = .true.  
    tstress       = .true.
    etot_conv_thr =  1.0d-6
    forc_conv_thr =  1.0d-5
   /

   &SYSTEM
    ibrav       = 0
    nat         = 2
    ntyp        = 1
    nbnd        = 16
    occupations = 'fixed'
    !    occupations = "smearing"
    !    smearing    = "cold"    
    !    degauss     =  1.0d-2
    ecutwfc     =  50
    ecutrho     =  400
    /

    &ELECTRONS
    conv_thr         =  1.000e-9
    electron_maxstep =  200
    mixing_beta      =  7.00000e-01
    startingpot      = "atomic"
    startingwfc      = "atomic+random"
    /

    &IONS
    ion_dynamics = "bfgs"
    /

    &CELL
    cell_dofree    = "x"
    cell_dynamics  = "bfgs"
    press_conv_thr =  0.01
    /

    K_POINTS {automatic}
     40  1  1  0 0 0

    ATOMIC_SPECIES
    C      12.01070  C.pbe-n-kjpaw_psl.1.0.0.UPF

    CELL_PARAMETERS (angstrom)
    2.565602620   0.000000000   0.000000000
    0.000000000  10.000000000   0.000000000
    0.000000000   0.000000000  10.000000000

    ATOMIC_POSITIONS (angstrom)
    C             0.0002913894       5.0000000000       5.0000000000
    C             1.2644833916       5.0000000000       5.0000000000
   ```

   采用下面命令查看计算过程中的受力情况已经压力张量和晶格常数变化情况  

   ```bash
   cat vc-relax.out |grep -A 10 "Total force ="
   ```

   计算结束后，使用下面命令得到优化后的结果

   ```bash
   awk  '/Begin final coordinates/,/End final coordinates/{print $0}' vc-relax.out

   
   Begin final coordinates
        new unit-cell volume =   1731.62982 a.u.^3 (   256.60106 Ang^3 )
        density =      0.15545 g/cm^3
   
   CELL_PARAMETERS (angstrom)
      2.566010647   0.000000000   0.000000000
      0.000000000  10.000000000   0.000000000
      0.000000000   0.000000000  10.000000000
   
   ATOMIC_POSITIONS (angstrom)
   C             0.0019333750        5.0000000000        5.0000000000
   C             1.2630425527        5.0000000000        5.0000000000
   End final coordinates
      
   ```

   2.2 `cp vc-relax.in relax.in` 将上面优化得到的`CELL_PARAMETERS` 和`ATOMIC_POSITIONS` 结果在relax.in中进行修改，并修改为`calculation = "relax"`, 将`&CELL /`部分注释掉。再对原子位置进行优化。计算完成后，运行下列命令得到优化后的原子位置。  

   `awk  '/Begin final coordinates/,/End final coordinates/{print $0}' relax.out`  
   得到如下结果：  

   ```bash
   Begin final coordinates

    ATOMIC_POSITIONS (angstrom)
    C             0.0019333750        5.0000000000        5.0000000000
    C             1.2630425527        5.0000000000        5.0000000000
    End final coordinates
   ```  

   2.3 `cp relax.in scf.in` 将上面的原子位置（**ATOMIC_POSITIONS**）结果在scf.in文件中进行修改。修改`calculation = "scf"` 将`&IONS /` 部分注释掉。`mkdir ../scf` 并将修改后的scf.in文件复制到`../scf`中进行自洽计算。  

   2.4 计算态密度.修改`calculaiton='nscf'` 并增加`kpoint`的采样密度。

   ```bash
   cp -r scf dos
   cd dos
   mv scf.in nscf.in
   ```  

   使用dos.x对nscf计算得态密度进行处理。输入文件`dos.in`如下：  

   ```fortran
   &DOS
   prefix = 'carbyne'
   outdir = './outdir'
   bz_sum = "smearing"
   ngauss = 0
   degauss = 2.0d-2
   DeltaE  = 0.01
   fildos = 'carbyne.dos'
   /
   ```

   对计算出来的态密度采用OriginPro进行作图，(采用较密kpoint计算的结果)如下所示。  
   ![dos](https://github.com/xh125/MarkdownImage/raw/main/Image/LVCSH/dos.png)  
   采用 **`projwfc.x`** 对态密度进行分波态密度计算，用于bandfat分析，后面进行Wannier计算，需要根据分波态密度来选取初始猜测的wannier函数。  

   ```fortran
   projwfc.in
    
    &PROJWFC
   prefix = 'carbyne'
   outdir = './outdir'
   ngauss = 0
   degauss= 0.01
   DeltaE = 0.5
   filpdos= 'carbyne.pdos'
   filproj= 'carbyne.proj'
   / 
   
   projwfc.x <projwfc.in>projwfc.out
   ```  

   处理分波态密度得到的结果如下:  

   ![Pdos](https://github.com/xh125/MarkdownImage/raw/main/Image/LVCSH/Pdos.png)  

   2.5 计算能带,修改`calculation='bands'`，`nbnd` 和`K_POINTS crystal_b`部分进行第一布里渊区高对称kpoint的能带计算。  

   ```bash
   cp -r scf band
   cd band
   cp scf.in pw-bands.in
   
   calculation   = "bands"
   nbnd = 22
   K_POINTS crystal_b
   3
   -0.5 0.0 0.0 50
    0.0 0.0 0.0 50
    0.5 0.0 0.0 1
   ```

   计算完成后使用bands.x处理能带数据。在能带计算之后，用projwfc.x生成的fatband.projwfc_up文件，可以和bands.x生成的文件bands.dat结合，画出各个能带的原子轨道投影，画图脚本如下：/LVCSH/tools/fatband.f90。计算能带图如下:  
   ![band.png](https://github.com/xh125/MarkdownImage/raw/main/Image/LVCSH/fatband.png)  

   2.6 根据fatband中的结果，进行wannier90计算  

   2.6.1 Scf计算  

   2.6.2 nscf计算

   ```bash
   cp scf.in nscf.in并做如下修改
   calculation   = "nscf" 
   nbnd        = 16
   !K_POINTS {automatic}
   ! 1  1  40  0 0 0
   
   kmesh.pl 1 1 40 >>nscf.in
   pw.x < nscf.in > nscf.out
   
   ```  

   2.6.3 Wannier90计算  

   ```fortran
   !System
   num_wann  = 4
   num_bands = 6
    
   !Projection
   begin projections
   C : px;py
   end projections
    
   !Job Control
   exclude_bands : 1-2
   !restart =
    
   !disentanglement
   dis_win_min  = -25.0
   dis_win_max  = 7.0
   dis_froz_min = -12.0
   dis_froz_max = -0.6
   dis_num_iter = 1000
   dis_mix_ratio= 0.5
   dis_conv_tol = 1.0E-10
   dis_conv_window = 5
    
   !Wannierise
   num_iter = 10000
   conv_tol = 1.0E-10
   conv_window     = 10
   guiding_centres = .true.

   ! SYSTEM
    
   begin unit_cell_cart
   Ang
     20.000000   0.000000   0.000000
     0.000000  20.000000   0.000000
     0.000000   0.000000   2.565985410
   end unit_cell_cart
    
   begin atoms_cart
   Ang
   C       10.0000000000       10.0000000000       -0.0019756118
   C       10.0000000000       10.0000000000        1.2517014809
   end atoms_cart
    
   ! KPOINTS     
   ```  

   ```bash
   echo "mp_grid = 1 1 40">>${seedname}.win
   echo "begin kpoints">>${seedname}.win
   kmesh.pl 1 1 40 wannier >>${seedname}.win
   echo "end kpoints">>${seedname}.win
   ```  

   2.6.4 'wannier90.x -pp seedname'  
   2.6.5 `pw2wannier90.x < pw2wan.in > pw2wan.out`  

   ```fortran
   pw2wan.in
   &inputpp
      outdir   = './outdir'
      prefix   = 'carbyne'
      seedname = 'carbyne'
      spin_component = 'none'
      write_mmn = .true.
      write_amn = .true.
      write_unk = .true.
   /
   ```  

   2.6.6 `mpirun -np $NP wannier90.x seedname`

>In directory epw to calculate the electron-phonon coupling matrix using the changed EPW code. And the output be named dependend on the kpoint: as epw40.out, epw80.out, epw120.out, epw160.out. Used to test the kpoint and qpoint convergence.  

1. make a directory for lvcsh calculation  

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

5. Use shell script mkepwdir.sh to build dir for diffrent kpoints and qpoints Surface hopping calculation. And the script will make a test running in the different director to give how to set **`nefre_sh`** and **`nhfre_sh`**.  

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
        sed -i "s/lvcsh-epw/lvcsh-epw$i-n0/g" epw$i/lvcsh.bsub
        cp job.sh epw$i
        cp LVCSH.in epw$i
        sed -i "s:epw40:epw$i:g" epw$i/LVCSH.in
        sed -i "s:ncore:ncore=$ncore #:g" epw$i/LVCSH.in
        cp LVCSH.in epw$i/QEfiles
        sed -i "s/low/high/g" epw$i/QEfiles/LVCSH.in
        sed -i "s:../../QEfiles/epw40.out:epw$i.out:g" epw$i/QEfiles/LVCSH.in
        cp lvcsh-test.bsub epw$i/QEfiles
        cd epw$i/QEfiles
        sed -i "2s/lvcsh-epw/lvcsh-epw$i/g" lvcsh-test.bsub
        bsub < lvcsh-test.bsub
        cd ../..
    done
```  

```bash
job.sh
#!/bin/bash
    nnodes=10
    for i in {1..$nnodes}
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
epwoutname    = "./QEfiles/epw40.out"
!nefre_sh      = 40
!nhfre_sh      = 40
nnode         = 1
ncore         = 28
naver         = 10
nsnap         = 1000
nstep         = 2
dt            = 0.5
savedsnap     = 1
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
```  

```bash
lvcsh-test.bsub
#!/bin/bash
#BSUB -J lvcsh-epw-test
#BSUB -q privateq-zw
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -o %J.out
#BSUB -e %J.err

export MODULEPATH=/share/home/zw/xiehua/opt/modules-4.7.1/modulefiles
module load lvcsh/0.6.2

LVCSH_complex.x

```  

6. By look the initial adiabatic state in the QEfiles/LVCSH.out for different kpoints directory. Set the **`nefre_sh`** and **`nhfre_sh`** in the QEfiles/LVCSH.in to tests the time for one step nonadiabatic calculation. Then, subscrib the job again.  

```shell
bsub < lvcsh-test.bsub
```

7. change the LVCSH.in file in the epw40, including the parameters for lvcsh run as following:

```fortran
nefre_sh = 34
nhfre_sh = 34
naver    = 100
nstep    = 2
nsnap    = 1000
dt       = 0.5
nnode    = 10
ncore    = 28
savedsnap= 25
```

change the **`nnodes`** in the job.sh bash script. Then './job.sh' to make lvcsh running in `nnodes` node, which include `ncore` CPU cores.  

8. After the setp 7, change the parameter `calculation = plot` and `epwoutname = './QEfiles/epw40.out'` of  LVCSH.in in the dir `epw40`. Then run `LVCSH_complex.x` in `epw40` dir. The code will read the result in all nodes and cores, then get a average results and writing in the dir `epw40`.

## Butterfly of the code

![Alt Butterfly](https://github.com/xh125/MarkdownImage/raw/main/Image/LVCSH/Butterfly-lvcsh.png)
