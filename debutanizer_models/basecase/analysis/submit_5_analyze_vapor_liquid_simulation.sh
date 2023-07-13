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

which python

PFM_PATH="/home/gridsan/hwpang/Software/PolymerFoulingModeling"
alpha_rates_path="simulation_results/1,3-BUTADIENE_1.0_3600.0_64.0/alpha_rates.yml"
simulation_results_path="simulation_results/1,3-BUTADIENE_1.0_3600.0_64.0/simulation_vapor_liquid_yliqn_3648.0.csv"
aspen_condition_path=$PFM_PATH/debutanizer_models/basecase/aspen_simulation/aspen_conditions.yml
model_name="basecase_debutanizer_model"

start=$(date +%s)
python-jl $PFM_PATH/debutanizer_models/basecase/analysis/analyze_vapor_liquid_simulation.py \
    --alpha_rates_path $alpha_rates_path \
    --model_name $model_name \
    --aspen_condition_path $aspen_condition_path \
    --simulation_results_path $simulation_results_path
end=$(date +%s)
runtime=$((end - start))
echo "runtime: $runtime"
