#!/bin/sh -l
#SBATCH -J analyze_film_simulation
#SBATCH -o slurm-analyze_film_simulation.out

conda activate rmg_py3_20230404

which julia
which python

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
model_name=trace_oxygen_perturbed_debutanizer_model
simulation_directory="../simulation_results/OXYGEN_1e0_3600.0_32.0"

start=$(date +%s)
python $PFM_PATH/debutanizer_models/trace_oxygen_perturbed/analysis/analyze_film_simulation.py \
        --model_name $model_name \
        --simulation_directory $simulation_directory
end=$(date +%s)
runtime=$((end - start))
echo "runtime: $runtime"
