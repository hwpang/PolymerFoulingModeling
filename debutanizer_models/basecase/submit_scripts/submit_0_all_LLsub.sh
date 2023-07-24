#!/bin/bash -l

# LLsub ./submit_0_all_LLsub.sh [1,1,48] -q spot-xeon-p8
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

chemkin_path=$PFM_PATH/debutanizer_models/basecase/liquid_mechanism/chem_annotated.inp
species_dict_path=$PFM_PATH/debutanizer_models/basecase/liquid_mechanism/species_dictionary.txt
model_name=basecase_debutanizer_model

start_time=$(date +%s)
python-jl $PFM_PATH/shared/liquid_mechanism_parameters_update/update_liquid_mechanism_parameters.py \
    --chemkin_path $chemkin_path \
    --species_dict_path $species_dict_path \
    --n_jobs 48 \
    --model_name $model_name
end_time=$(date +%s)

echo execution time was $(expr $end_time - $start_time) s.

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
chemkin_path=liquid_mechanism/chem_annotated_updated.inp
species_dict_path=liquid_mechanism/species_dictionary_updated.txt
model_name="basecase_debutanizer_model"

python-jl $PFM_PATH/shared/film_mechanism_generation/generate_film_mechanism.py \
    --chemkin_path $chemkin_path \
    --species_dict_path $species_dict_path \
    --n_jobs 48 \
    --model_name $model_name

conda deactivate

echo "============================================================"
echo "submit liquid simulation"
LLsub ./submit_0.1_liquid_simulation_LLsub.sh [18,1,48] -q spot-xeon-p8
echo "============================================================"
