#!/bin/bash -l

# LLsub ./submit_6_generate_film_mechanism_LLsub.sh [1,1,48] -q spot-xeon-p8
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
which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
chemkin_path=liquid_mechanism/chem_annotated.inp
species_dict_path=liquid_mechanism/species_dictionary.txt
model_name="trace_oxygen_perturbed_debutanizer_model"

python-jl $PFM_PATH/shared/film_mechanism_generation/generate_film_mechanism.py \
    --chemkin_path $chemkin_path \
    --species_dict_path $species_dict_path \
    --n_jobs 48 \
    --model_name $model_name

conda deactivate
