#!/bin/sh -l
#SBATCH -J analyze_sensitivity
#SBATCH -o slurm-analyze_sensitivity.out

conda activate rmg_py3_20230404

which julia
which python

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
model_name=trace_oxygen_perturbed_debutanizer_model
all_simulation_directory="../simulation_results"

start=$(date +%s)
python $PFM_PATH/debutanizer_models/trace_oxygen_perturbed/analysis/analyze_film_simulation.py \
        --model_name $model_name \
        --all_simulation_directory $all_simulation_directory
end=$(date +%s)
runtime=$((end - start))
echo "runtime: $runtime"
