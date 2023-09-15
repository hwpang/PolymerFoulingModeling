#!/bin/bash -l

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

conda activate rmg_py3_20230404

which julia

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
rms_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/liquid_mechanism/chem.rms
model_name="trace_oxygen_perturbed_debutanizer_model"

julia $PFM_PATH/shared/analysis/analyze_debutanizer_time_step_studies.jl \
        $rms_path \
        $model_name
        