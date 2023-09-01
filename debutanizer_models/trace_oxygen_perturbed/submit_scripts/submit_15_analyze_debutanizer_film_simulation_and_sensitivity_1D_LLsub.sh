#!/bin/bash -l

# LLsub ./submit_13_analyze_debutanizer_film_simulation_and_sensitivity_1D_LLsub.sh [1,1,48] -q spot-xeon-p8
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

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
model_name="trace_oxygen_perturbed_debutanizer_model"
simulation_directory="simulation_results/OXYGEN_1e0_3600.0_32.0"

start=$(date +%s)
python $PFM_PATH/shared/analysis/analyze_debutanizer_film_simulation_1D.py \
    --model_name $model_name \
    --simulation_directory $simulation_directory
end=$(date +%s)
runtime=$((end - start))
echo "runtime: $runtime"

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

conda activate rmg_py3_20230404

which python

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
model_name=trace_oxygen_perturbed_debutanizer_model
all_simulation_directory="simulation_results"

start=$(date +%s)
python $PFM_PATH/shared/analysis/analyze_debutanizer_sensitivity_1D.py \
    --model_name $model_name \
    --all_simulation_directory $all_simulation_directory
end=$(date +%s)
runtime=$((end - start))
echo "runtime: $runtime"