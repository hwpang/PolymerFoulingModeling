#!/bin/bash -l
#SBATCH -J generate_film_mechanism
#SBATCH -o slurm-generate_film_mechanism.out
#SBATCH -n 48

conda activate rmg_py3_20230404

which julia
which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
chemkin_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/liquid_mechanism/chem_annotated.inp
species_dict_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/liquid_mechanism/species_dictionary.txt
model_name="trace_oxygen_perturbed_debutanizer_model"

python-jl $PFM_PATH/debutanizer_models/basecase/film_mechanism_generation/generate_film_mechanism.py \
            --chemkin_path $chemkin_path \
            --species_dict_path $species_dict_path \
            --n_jobs 48 \
            --model_name $model_name

conda deactivate
