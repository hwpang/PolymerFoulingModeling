#!/bin/bash

# LLsub ./submit_3_generate_liquid_mechanism_LLsub.sh [1,12,4] -q spot-xeon-p8
# squeue -p spot-xeon-p8
# watch LLloadSpot

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

# Run job

conda activate rmg_py3_20230404
which julia
which python
which python-jl

RMG_PATH=/home/gridsan/hwpang/Software/RMG-Py
export PYTHONPATH=$RMG_PATH:$PYTHONPATH
export PATH=$RMG_PATH:$PATH

python-jl $RMG_PATH/rmg.py -n 12 input.py

conda deactivate




