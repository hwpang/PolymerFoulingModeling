#!/bin/bash -l

# LLsub ./submit_14_plotter.sh [1,1,48] -q spot-xeon-p8
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

# start=$(date +%s)
# python $PFM_PATH/shared/analysis/plot_rates.py
# end=$(date +%s)
# runtime=$((end - start))
# echo "runtime: $runtime"

start_time=$(date +%s)
python $PFM_PATH/shared/analysis/plot_thermo_GAV_comparison.py \
    --model_name trace_oxygen_perturbed_debutanizer_model
end_time=$(date +%s)
echo execution time was $(expr $end_time - $start_time) s.

start_time=$(date +%s)
python $PFM_PATH/shared/analysis/plot_model_size_and_parameter_sources.py \
    --model_name $model_name
end_time=$(date +%s)
echo execution time was $(expr $end_time - $start_time) s.
