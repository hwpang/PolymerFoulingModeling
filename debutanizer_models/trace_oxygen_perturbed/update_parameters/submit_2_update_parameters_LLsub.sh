#!/bin/bash -l

# LLsub ./submit_2_update_parameters_LLsub.sh [1,1,48] -q spot-xeon-p8
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
which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling

chemkin_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/liquid_mechanism/chem_annotated.inp
species_dict_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/liquid_mechanism/species_dictionary.txt

start_time=$(date +%s)
python-jl $PFM_PATH/utils/update_parameters/update_parameters.py \
    --chemkin_path $chemkin_path \
    --species_dict_path $species_dict_path \
    --n_jobs 48
end_time=$(date +%s)

echo execution time was $(expr $end_time - $start_time) s.
