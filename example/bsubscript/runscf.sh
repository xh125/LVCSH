#!/bin/bash
mkdir scf
cp -r pseudo scf/
cd scf
press=100
a1=$(grep -B 30 "Begin final coordinates" ../relax/vc-relax.out |grep "celldm(1)"|awk '{print $3}')
a3=$(grep -B 30 "Begin final coordinates" ../relax/vc-relax.out |grep "celldm(3)"|awk '{print $3}')
position=$(awk  '/Begin final coordinates/,/End final coordinates/{print $0}' ../relax/vc-relax.out|grep -A 8 "ATOMIC_POSITIONS (crystal)")
cat > scf.in << EOF
&CONTROL 
calculation   = "scf" 
prefix    = zngep2 
outdir    = './' 
pseudo_dir    = './pseudo' 
verbosity = 'high' 
tprnfor   = .true.
tstress   = .true.
etot_conv_thr =  1.00000e-05 
forc_conv_thr =  1.00000e-04 
nstep     = 100 
!    restart_mode  = "from_scratch" 
/  
&SYSTEM 
ibrav   =   7 
celldm(1)   =   $a1
celldm(3)   =   $a3   
nat     =   8 
ntyp    =   3 
ecutwfc =   90 
!    ecutrho =   300  
occupations = 'fixed' 
!    occupations = "smearing" 
!    smearing    = "gaussian" 
!    degauss = 0.01 
/  
&ELECTRONS 
conv_thr     =  1.00000e-08 
electron_maxstep =  200 
mixing_beta  =  0.7 
diagonalization  =  'david' 
!    startingpot  = "atomic" 
!    startingwfc  = "atomic+random" 
/  
!&IONS 
!    ion_dynamics = "bfgs" 
!/  
&CELL 
!    cell_dofree    = "ibrav" 
!    cell_dynamics  = "bfgs" 
press=$press 
!    press_conv_thr =  0.01 
/  
K_POINTS {automatic} 
10 10 10  0 0 0  
ATOMIC_SPECIES 
Zn 65.39000  Zn_ONCV_PBE-1.2.upf 
Ge 72.61000  Ge_ONCV_PBE-1.2.upf 
P  30.97376  P_ONCV_PBE-1.2.upf 

$position
EOF

cat > qe-scf.bsub << EOF
#!/bin/bash 
#BSUB -J qe-scf 
#BSUB -q publicq
#BSUB -n 28 
#BSUB -R "span[ptile=28]" 
#BSUB -o %J.out 
#BSUB -e %J.err  

CURDIR=\$PWD 
#Generate nodelist 
rm -f \${CURDIR}/nodelist >& /dev/null 
for i in `echo \$LSB_HOSTS`; do 
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
source /share/home/zw/xiehua/opt/modules-5.0.0/init/profile.sh export 
MODULEPATH=/share/home/zw/xiehua/modulefiles  

module load Quantum_Espresso/6.8.0

mpirun -np \$NP pw.x -nk 7 <scf.in> scf.out 

EOF

bsub < qe-scf.bsub
