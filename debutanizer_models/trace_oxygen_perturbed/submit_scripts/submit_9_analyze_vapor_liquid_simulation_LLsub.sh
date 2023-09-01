#!/bin/bash -l

# LLsub ./submit_5_analyze_vapor_liquid_simulation_LLsub.sh [1,1,48] -q spot-xeon-p8
# squeue -p spot-xeon-p8
# watch LLloadSpot

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
time_step=32.0
simulation_results_directory="simulation_results/OXYGEN_1e0_3600.0_${time_step}"
aspen_condition_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/aspen_simulation/aspen_conditions.yml
model_name="trace_oxygen_perturbed_debutanizer_model"

start=$(date +%s)
python $PFM_PATH/shared/analysis/analyze_debutanizer_vapor_liquid_simulation.py \
    --model_name $model_name \
    --aspen_condition_path $aspen_condition_path \
    --simulation_results_directory $simulation_results_directory
end=$(date +%s)
runtime=$((end - start))
echo "runtime: $runtime"
