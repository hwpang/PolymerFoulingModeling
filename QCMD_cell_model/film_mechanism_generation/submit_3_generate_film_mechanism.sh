#!/bin/bash -l
#SBATCH -o slurm-%x.out
#SBATCH -n 48
#SBATCH -N 1

conda activate rmg_py3_20230404

which julia
which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
chemkin_path=$PFM_PATH/QCMD_cell_model/liquid_mechanism/chem_annotated.inp
species_dict_path=$PFM_PATH/QCMD_cell_model/liquid_mechanism/species_dictionary.txt
model_name="QCMD_cell_model"

python-jl $PFM_PATH/debutanizer_models/basecase/film_mechanism_generation/generate_film_mechanism.py \
    --chemkin_path $chemkin_path \
    --species_dict_path $species_dict_path \
    --n_jobs 48 \
    --model_name $model_name

conda deactivate
