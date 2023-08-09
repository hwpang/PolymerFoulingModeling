#!/bin/bash -l

# LLsub ./submit_1_extract_and_plot_aspen_conditions.sh [1,1,48] -q spot-xeon-p8
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

which julia
which python

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling

model_name=trace_oxygen_perturbed_debutanizer_model
aspen_results_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/aspen_simulation/Distillation_FEDP_XN_V2_060720_v2_oxygen.xlsx

start_time=$(date +%s)
python $PFM_PATH/debutanizer_models/basecase/aspen_simulation/extract_and_plot_aspen_conditions.py
--model_name $model_name \
    --aspen_results_path $aspen_results_path
end_time=$(date +%s)

echo execution time was $(expr $end_time - $start_time) s.
