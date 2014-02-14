#!/bin/bash -l
#SBATCH --job-name="multigpu_lbm" 
#SBATCH --nodes=256
#SBATCH --ntasks=256
#SBATCH --cpus-per-task=1 
#SBATCH --ntasks-per-node=1
#SBATCH --time=11:00:00 
#======START-BASIC=============================== 
NGPUs=256
NEXPs=2
# module swap PrgEnv-cray PrgEnv-intel
# git checkout basic
# make clean
# scons --compiler=CC benchmark=1 --mode=release
./benchmark.py $NGPUs $NEXPs aprun ./build/lbm_opencl_dc_CC_release weak
mkdir -p ./output/results_benchmark/todi/weak/1d/basic/
mkdir -p ./output/benchmark/
mv -f ./output/benchmark/*.ini ./output/results_benchmark/todi/weak/1d/basic/
mv -f ./output/benchmark/*.txt ./output/results_benchmark/todi/weak/1d/basic/
mv -f ./output/benchmark/*.out ./output/results_benchmark/todi/weak/1d/basic/
#======END=================================
# #======START-=============================== 
# git checkout ccol_sbk_scq
# make clean
# scons --compiler=CC benchmark=1 --mode=release
# ./benchmark.py $NGPUs $NEXPs aprun ./build/lbm_opencl_dc_CC_release weak
# mkdir -p ./output/results_benchmark/todi/weak/1d/ccol_sbk_scq
# mv -f ./output/benchmark/*.ini ./output/results_benchmark/todi/weak/1d/ccol_sbk_scq
# mv -f ./output/benchmark/*.txt ./output/results_benchmark/todi/weak/1d/ccol_sbk_scq
# mv -f ./output/benchmark/*.out ./output/results_benchmark/todi/weak/1d/ccol_sbk_scq
# #======END=================================
# #======START-=============================== 
# git checkout ccol_mbk_scq_nocb
# make clean
# scons --compiler=CC benchmark=1 --mode=release
# ./benchmark.py $NGPUs $NEXPs aprun ./build/lbm_opencl_dc_CC_release weak
# mkdir -p ./output/results_benchmark/todi/weak/1d/ccol_mbk_scq_nocb
# mv -f ./output/benchmark/*.ini ./output/results_benchmark/todi/weak/1d/ccol_mbk_scq_nocb
# mv -f ./output/benchmark/*.txt ./output/results_benchmark/todi/weak/1d/ccol_mbk_scq_nocb
# mv -f ./output/benchmark/*.out ./output/results_benchmark/todi/weak/1d/ccol_mbk_scq_nocb
# #======END=================================
# #======START-=============================== 
# git checkout ccol_mbk_mcq_nocb
# make clean
# scons --compiler=CC benchmark=1 --mode=release
# ./benchmark.py $NGPUs $NEXPs aprun ./build/lbm_opencl_dc_CC_release weak
# mkdir -p ./output/results_benchmark/todi/weak/1d/ccol_mbk_mcq_nocb
# mv -f ./output/benchmark/*.ini ./output/results_benchmark/todi/weak/1d/ccol_mbk_mcq_nocb
# mv -f ./output/benchmark/*.txt ./output/results_benchmark/todi/weak/1d/ccol_mbk_mcq_nocb
# mv -f ./output/benchmark/*.out ./output/results_benchmark/todi/weak/1d/ccol_mbk_mcq_nocb
# #======END=================================
