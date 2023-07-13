#!/bin/bash -l
#SBATCH -n 48
#SBATCH -o slurm-%x.out
#SBATCH -N 1

conda activate rmg_py3_20230404

which julia
which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling

chemkin_path=$PFM_PATH/QCMD_cell_model/liquid_mechanism/chem_annotated.inp
species_dict_path=$PFM_PATH/QCMD_cell_model/liquid_mechanism/species_dictionary.txt

python-jl $PFM_PATH/utils/update_parameters/update_parameters.py \
    --chemkin_path $chemkin_path \
    --species_dict_path $species_dict_path \
    --n_jobs 48
