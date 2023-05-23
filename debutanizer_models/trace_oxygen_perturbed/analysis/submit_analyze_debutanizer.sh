#!/bin/bash -l
#SBATCH -J analyze_debutanizer
#SBATCH -o slurm-analyze_debutanizer.out

conda activate rmg_py3_20230404

which python

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
alpha_rates_path="../simulation_results/OXYGEN_1e0_3600.0_32.0/alpha_rates.yml"
results_path="../simulation_results/OXYGEN_1e0_3600.0_32.0/simulation_vapor_liquid_yliqn_3616.0.csv"
model_name=trace_oxygen_perturbed_debutanizer_model

start=$(date +%s)
python $PFM_PATH/debutanizer_models/trace_oxygen_perturbed/analysis/analyze_debutanizer_model.py \
        --alpha_rates_path $alpha_rates_path \
        --model_name $model_name \
        --results_path $results_path
end=$(date +%s)
runtime=$((end - start))
echo "runtime: $runtime"
