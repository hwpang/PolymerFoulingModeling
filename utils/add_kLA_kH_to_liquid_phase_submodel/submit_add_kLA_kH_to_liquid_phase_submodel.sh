#!/bin/bash -l
#SBATCH -J add_kLA_kH_to_liquid_phase_submodel
#SBATCH -o slurm-add_kLA_kH_to_liquid_phase_submodel.out

conda activate rmg_py3_20230404

which julia
which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling

chemkin_path=$1
species_dict_path=$2
rms_path=$3

python-jl $PFM_PATH/utils/add_kLA_kH_to_liquid_phase_submodel/add_kLA_kH_to_liquid_phase_submodel.py \
            --chemkin_path $chemkin_path \
            --species_dict_path $species_dict_path \
            --rms_path $rms_path