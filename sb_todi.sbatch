#!/bin/bash -l
#SBATCH --job-name="multigpu_lbm" 
#SBATCH --nodes=272
#SBATCH --ntasks=272
#SBATCH --cpus-per-task=1 
#SBATCH --ntasks-per-node=1
#SBATCH --time=03:00:00 
#======START=============================== 
echo "On which nodes it executes" 
echo $SLURM_JOB_NODELIST 
echo "Now run the MPI tasks..." 
./benchmark.py 272 1 aprun ./build/lbm_opencl_dc_CC_release weak
#======END=================================
