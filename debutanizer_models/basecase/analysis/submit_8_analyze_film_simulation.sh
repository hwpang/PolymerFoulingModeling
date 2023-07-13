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

which julia
which python

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
model_name=basecase_debutanizer_model
simulation_directory="simulation_results/1,3-BUTADIENE_1.0_3600.0_64.0"

start=$(date +%s)
python $PFM_PATH/debutanizer_models/basecase/analysis/analyze_film_simulation.py \
    --model_name $model_name \
    --simulation_directory $simulation_directory
end=$(date +%s)
runtime=$((end - start))
echo "runtime: $runtime"
