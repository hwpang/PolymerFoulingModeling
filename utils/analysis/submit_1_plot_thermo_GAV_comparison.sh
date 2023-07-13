#!/bin/bash -l
#SBATCH -o slurm-%x.out

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

conda activate rmg_py3_20230404

which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling

start_time=$(date +%s)
python-jl $PFM_PATH/utils/analysis/plot_thermo_GAV_comparison.py
end_time=$(date +%s)

echo execution time was $(expr $end_time - $start_time) s.
