#!/bin/bash -l
#SBATCH -o slurm-%x-%a.out
#SBATCH -c 2
#SBATCH --array=0-150

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

conda activate rmg_py3_20230404

which julia
which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling

rms_mech_directory=$PFM_PATH/debutanizer_models/basecase/film_mechanism
aspen_condition_path=$PFM_PATH/debutanizer_models/basecase/aspen_simulation/aspen_conditions.yml
model_name="basecase_debutanizer_model"

jobs=()
perturb_species_list=("1,3-BUTADIENE" "CYCLOPENTADIENE")
perturb_factor_list=("0.5" "0.7" "0.9" "1.0" "1.1" "1.3" "1.5" "1.7" "1.9")

for perturb_species in "${perturb_species_list[@]}"; do
    for perturb_factor in "${perturb_factor_list[@]}"; do
        for tray in {1..40}; do
            jobs+=("$perturb_species $perturb_factor $tray")
        done
    done
done

for jobind in $(seq $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT ${#jobs[@]}); do
    echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
    echo "SLURM_ARRAY_TASK_COUNT: $SLURM_ARRAY_TASK_COUNT"
    echo "jobind: $jobind"
    echo "${jobs[$jobind]}"
    job=(${jobs[$jobind]})
    perturb_species=${job[0]}
    perturb_factor=${job[1]}
    tray=${job[2]}
    start=$(date +%s)
    liquid_simulation_results_path="simulation_results/${perturb_species}_${perturb_factor}_3600.0_64.0/simulation_vapor_liquid_yliqn_3648.0.csv"
    julia $PFM_PATH/debutanizer_models/basecase/film_simulation/simulate_film_submodel.jl \
        $rms_mech_directory \
        $model_name \
        $liquid_simulation_results_path \
        $aspen_condition_path \
        $tray
    end=$(date +%s)
    runtime=$((end - start))
    echo "runtime: $runtime"
done
