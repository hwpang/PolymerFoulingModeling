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
model_name=basecase_debutanizer_model
all_simulation_directory="simulation_results"

start=$(date +%s)
python-jl $PFM_PATH/debutanizer_models/basecase/analysis/analyze_sensitivity.py \
    --model_name $model_name \
    --all_simulation_directory $all_simulation_directory
end=$(date +%s)
runtime=$((end - start))
echo "runtime: $runtime"
