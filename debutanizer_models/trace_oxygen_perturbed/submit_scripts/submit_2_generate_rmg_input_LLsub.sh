#!/bin/bash -l

# LLsub ./submit_2_generate_rmg_input_LLsub.sh [1,1,48] -q spot-xeon-p8
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
which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
aspen_condition_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/aspen_simulation/aspen_conditions.yml
python-jl $PFM_PATH/shared/liquid_mechanism_rmg_input_generation/generate_liquid_mechanism_rmg_input.py \
            --aspen_condition_path $aspen_condition_path \
            --model_name trace_oxygen_perturbed_debutanizer_model