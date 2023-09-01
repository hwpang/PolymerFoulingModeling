#!/bin/bash -l

# LLsub ./submit_4_get_alpha_for_ASF_distribution_LLsub.sh [8,1,48] -q spot-xeon-p8
# squeue -p spot-xeon-p8
# watch LLloadSpot

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

chemkin_path=liquid_mechanism/chem_annotated.inp
species_dict_path=liquid_mechanism/species_dictionary.txt
rms_path=liquid_mechanism/chem.rms
aspen_condition_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/aspen_simulation/aspen_conditions.yml
model_name="trace_oxygen_perturbed_debutanizer_model"

jobs=()
perturb_species_list=("OXYGEN")
perturb_factor_list=("0.0" "1e-6" "1e-5" "1e-4" "1e-3" "1e-2" "1e-1" "1e0")

for perturb_species in "${perturb_species_list[@]}"; do
    for perturb_factor in "${perturb_factor_list[@]}"; do
        jobs+=("$perturb_species $perturb_factor")
    done
done

for jobind in $(seq $LLSUB_RANK $LLSUB_SIZE ${#jobs[@]}); do
    echo "LLSUB_RANK: $LLSUB_RANK"
    echo "LLSUB_SIZE: $LLSUB_SIZE"
    echo "jobind: $jobind"
    echo "${jobs[$jobind]}"
    job=(${jobs[$jobind]})
    perturb_species=${job[0]}
    perturb_factor=${job[1]}
    start=$(date +%s)
    time_step="32.0"
    liquid_simulation_results_directory="simulation_results/${perturb_species}_${perturb_factor}_3600.0_${time_step}"
    python-jl $PFM_PATH/debutanizer_models/basecase/ASF_distribution/get_alpha_for_ASF_distribution.py \
        --chemkin_path $chemkin_path \
        --species_dict_path $species_dict_path \
        --rms_path $rms_path \
        --aspen_condition_path $aspen_condition_path \
        --results_directory $liquid_simulation_results_directory \
        --n_jobs 48 \
        --model_name $model_name \
        --time_step $time_step
    end=$(date +%s)
    runtime=$((end - start))
    echo "runtime: $runtime"
done
