#!/bin/sh
#$ -N stars_only_run 
#$ -M hr203@sussex.ac.uk 
#$ -m bea
#$ -cwd
#$ -pe openmpi 104#256 #128 


#$ -q eng-inf_parallel.q #- not allowed to use this anymore
##$ -q parallel.q 
## #$ -q mps.q
## this hardware is not really suitable for MPI
##$ -q mps.q@@compute_amd 

##$ -q mps.q@@compute_intel
##$  -q inf_amd.q@@compute_amd 


#$ -S /bin/bash
# source modules environment:
module add sge #for submitting jobs
module load intel
module load mvapich2/intel/64/1.9 

#module load gcc
#module load openmpi/gcc/64/1.7.3
export OMP_NUM_THREADS=8
#which mpirun
#mpirun -np $NSLOTS ./C2Ray_3D_test_kyl_periodic_omp testinput
mpirun -np $NSLOTS ./c_alpha

