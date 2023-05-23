#!/bin/sh -l
#SBATCH -J analyze_film_simulation
#SBATCH -o slurm-analyze_film_simulation.out

conda activate rmg_py3_20230404

which julia

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
model_name=trace_oxygen_perturbed_debutanizer_model
rms_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/liquid_mechanism/chem.rms

start=$(date +%s)
julia $PFM_PATH/debutanizer_models/basecase/analysis/analyze_film_simulation.jl \
        $model_name \
        $rms_path
runtime=$((end - start))
echo "runtime: $runtime"