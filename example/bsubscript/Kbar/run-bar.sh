#!/bin/bash
for i in $(seq 300 20 300);do 
mkdir $i.Kbar

cd $i.Kbar
mkdir relax
cd relax
cp -r ../../pseudo ./

cat > vc-relax.in << EOF
&CONTROL
    calculation   = "vc-relax"
    prefix        = zngep2
    outdir        = './'
    pseudo_dir    = './pseudo'
    verbosity     = 'high'
    tprnfor       = .true.
    tstress       = .true.
    etot_conv_thr =  1.00000e-05
    forc_conv_thr =  1.00000e-04
    nstep         = 100
!    restart_mode  = "from_scratch"
/

&SYSTEM
    ibrav       =  7
    celldm(1) =     10.37600462
    celldm(3) =      1.97244178
    nat         =   8
    ntyp        =   3
    ecutwfc     =   90
!    ecutrho     =   300
     occupations = 'fixed'
!    occupations = "smearing"
!    smearing    = "gaussian"
!    degauss     = 0.01
/

&ELECTRONS
    conv_thr         =  1.00000e-08
    electron_maxstep =  200
    mixing_beta      =  0.7
    diagonalization  =  'david'
!    startingpot      = "atomic"
!    startingwfc      = "atomic+random"
/

&IONS
    ion_dynamics = "bfgs"
/

&CELL
    cell_dofree    = "ibrav"
    cell_dynamics  = "bfgs"
    press=$i
    press_conv_thr =  0.01
/

K_POINTS {automatic}
10 10 10  0 0 0

ATOMIC_SPECIES
Zn     65.39000  Zn_ONCV_PBE-1.2.upf
Ge     72.61000  Ge_ONCV_PBE-1.2.upf
P      30.97376  P_ONCV_PBE-1.2.upf

ATOMIC_POSITIONS (crystal)
Zn           -0.0000000000        0.0000000000        0.0000000000
Zn           -0.5000000000        0.7500000000        0.2500000000
Ge            0.5000000000        0.2500000000       -0.2500000000
Ge            0.0000000000        0.5000000000       -0.5000000000
P             0.0029929585        0.3750000000       -0.1279929585
P            -0.0029929585        0.8750000000       -0.6220070415
P             0.5029929585        0.6220070415       -0.3750000000
P            -0.5029929585        1.1279929585        0.1250000000

EOF

cat > qe-relax.bsub << EOF
#!/bin/bash
#BSUB -J $i.Kbar
#BSUB -q privateq-zw
#BSUB -n 28
#BSUB -R "span[ptile=28]"
#BSUB -o %J.out
#BSUB -e %J.err

CURDIR=\$PWD
#Generate nodelist
rm -f \${CURDIR}/nodelist >& /dev/null
for i in ``echo ` `\$LSB_HOSTS``
do
echo \$i >> \${CURDIR}/nodelist
done

NP=``cat ` `\${CURDIR}/nodelist` `|wc -l``
rm -f \${CURDIR}/nodelist >& /dev/null

#export OMP_NUM_THREADS=1
#export MKL_NUM_THREADS=1

# wgs 1
#source /share/home/ZhuangW/xh/opt/modules/init/profile.sh 
#export MODULEPATH=/share/home/ZhuangW/xh/modulefiles
# wgs 2
source /share/home/zw/xiehua/opt/modules-5.0.0/init/profile.sh
export MODULEPATH=/share/home/zw/xiehua/modulefiles

module load Quantum_Espresso/6.8.0

mpirun -np \$NP pw.x -nk 7 <vc-relax.in> vc-relax.out

EOF

bsub < qe-relax.bsub
cd ../..

done
