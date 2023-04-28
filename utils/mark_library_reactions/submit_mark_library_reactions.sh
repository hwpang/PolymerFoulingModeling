#!/bin/bash -l
#SBATCH -J mark_library_reactions
#SBATCH -o slurm-mark_library_reactions.out
#SBATCH -n 48

conda activate rmg_py3_20230404

which julia
which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling

chemkin_path=$1
species_dict_path=$2
rms_path=$3
model_name=$4

python-jl $PFM_PATH/utils/mark_library_reactions/mark_library_reactions.py \
            --chemkin_path $chemkin_path \
            --species_dict_path $species_dict_path \
            --rms_path $rms_path \
            --model_name $model_name \
            --n_jobs 48