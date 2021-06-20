#!/bin/bash
#PBS -N pw-vc-relax
#PBS -o ${PBS_JOBID}.out
#PBS -e ${PBS_JOBID}.err
##PBS -q fjwgs share
#PBS -q share
##PBS -l nodes=1:ppn=10
#PBS -l nodes=c3622:ppn=28

# get the number of processors
#NP=`cat $PBS_NODEFILE | wc -l`

# enter job's working directory
#cd $PBS_O_WORKDIR

#export OMP_NUM_THREADS=1
#export MKL_NUM_THREADS=1
export MODULEPATH=/public/software/modules:/public/home/wzhuang/xh/modulefiles

module load apps/Materials/Quantum_Espresso/6.7.0

export PBS_JOBNAME=pw-relax
export NP=14

# pw 第一个"-"之前的字符串
export exename1=${PBS_JOBNAME%%-*}.x
# scf 第一个“-”之后的所有字符串
export infilename1=${PBS_JOBNAME#*-}

#export exename2=ph.x
#export infilename2=ph

mpirun -np $NP ${exename1} -npool $NP <${infilename1}.in> ${infilename1}.out&
#mpirun -np $NP ${exename1} <${infilename1}.in> ${infilename1}.out
#mpirun -np $NP ${exename2} -npool $NP <${infilename2}.in> ${infilename2}.out

