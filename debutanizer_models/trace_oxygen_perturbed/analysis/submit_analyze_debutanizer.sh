#!/bin/bash -l
#SBATCH -J analyze_debutanizer
#SBATCH -o slurm-analyze_debutanizer.out

conda activate rmg_py3_20230404

which julia

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
alpha_rates_path="../simulation_results/OXYGEN_1e0_3600.0_32.0/alpha_rates.yml"
model_name=trace_oxygen_perturbed_debutanizer_model
julia $PFM_PATH/analyze_debutanizer_model.py \
        --alpha_rates_path $alpha_rates_path \
        --model_name $model_name
