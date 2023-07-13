#!/bin/bash -l
#SBATCH -J simulate_film_submodel
#SBATCH -o slurm-simulate_film_submodel-%a.out
#SBATCH --array=0-4

conda activate rmg_py3_20230404

which julia
which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling

rms_path=$PFM_PATH/QCMD_cell_model/film_mechanism/chem_film_phase.rms
model_name="QCMD_cell_model"

jobs=()
perturb_factor_list=("0.0" "1e-3" "1e-2" "1e-1" "1e0")

for perturb_factor in "${perturb_factor_list[@]}"; do
    jobs+=("$perturb_factor")
done

for jobind in $(seq $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT ${#jobs[@]}); do
    echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
    echo "SLURM_ARRAY_TASK_COUNT: $SLURM_ARRAY_TASK_COUNT"
    echo "jobind: $jobind"
    job=(${jobs[$jobind]})
    perturb_factor=${job[0]}
    echo "perturb_factor: $perturb_factor"
    liquid_simulation_results_path="simulation_results/O2_$perturb_factor/simulation_liquid_1.csv"
    echo "liquid_simulation_results_path: $liquid_simulation_results_path"
    start=$(date +%s)
    julia $PFM_PATH/debutanizer_models/basecase/film_simulation/simulate_film_submodel.jl \
        $rms_path \
        $model_name \
        $liquid_simulation_results_path
    end=$(date +%s)
    runtime=$((end - start))
    echo "runtime: $runtime"
done
