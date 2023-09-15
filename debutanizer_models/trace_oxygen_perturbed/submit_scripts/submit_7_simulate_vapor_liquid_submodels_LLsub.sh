#!/bin/bash -l

# LLsub ./submit_3_simulate_vapor_liquid_submodels_LLsub.sh [2,6,8] -q spot-xeon-p8
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

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling

aspen_condition_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/aspen_simulation/aspen_conditions.yml
rms_path=liquid_mechanism/chem_updated.rms
model_name=trace_oxygen_perturbed_debutanizer_model

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
    tray=${job[2]}
    start=$(date +%s)
    julia --threads 8 \
        $PFM_PATH/shared/debutanizer_vapor_liquid_simulation/simulate_vapor_liquid_submodels.jl \
        $aspen_condition_path \
        $rms_path \
        $perturb_species \
        $perturb_factor \
        $tray \
        32.0 \
        $model_name
    end=$(date +%s)
    runtime=$((end - start))
    echo "runtime: $runtime"
done
