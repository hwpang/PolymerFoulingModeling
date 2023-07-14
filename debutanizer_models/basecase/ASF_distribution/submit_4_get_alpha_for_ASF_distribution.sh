#!/bin/bash -l
#SBATCH -o slurm-%x-%a.out
#SBATCH -n 48
#SBATCH -N 1
#SBATCH --array=0-7

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

chemkin_path=$PFM_PATH/debutanizer_models/basecase/liquid_mechanism/chem_annotated.inp
species_dict_path=$PFM_PATH/debutanizer_models/basecase/liquid_mechanism/species_dictionary.txt
rms_path=$PFM_PATH/debutanizer_models/basecase/liquid_mechanism/chem.rms
aspen_condition_path=$PFM_PATH/debutanizer_models/basecase/aspen_simulation/aspen_conditions.yml
model_name="basecase_debutanizer_model"

jobs=()
perturb_species_list=("1,3-BUTADIENE" "CYCLOPENTADIENE")
perturb_factor_list=("0.5" "0.7" "0.9" "1.0" "1.1" "1.3" "1.5" "1.7" "1.9")

for perturb_species in "${perturb_species_list[@]}"; do
    for perturb_factor in "${perturb_factor_list[@]}"; do
        jobs+=("$perturb_species $perturb_factor")
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
    start=$(date +%s)
    liquid_simulation_results_directory="./simulation_results/${perturb_species}_${perturb_factor}_3600.0_64.0"
    python-jl $PFM_PATH/debutanizer_models/basecase/ASF_distribution/get_alpha_for_ASF_distribution.py \
        --chemkin_path $chemkin_path \
        --species_dict_path $species_dict_path \
        --rms_path $rms_path \
        --aspen_condition_path $aspen_condition_path \
        --results_directory $liquid_simulation_results_directory \
        --n_jobs 48 \
        --model_name $model_name
    end=$(date +%s)
    runtime=$((end - start))
    echo "runtime: $runtime"
done