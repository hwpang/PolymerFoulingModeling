#!/bin/bash -l
#SBATCH -o slurm-%x-%a.out
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --array=0-17

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

aspen_condition_path=$PFM_PATH/debutanizer_models/basecase/aspen_simulation/aspen_conditions.yml
rms_path=$PFM_PATH/debutanizer_models/basecase/liquid_mechanism/chem.rms

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
    tray=${job[2]}
    start=$(date +%s)
    julia --threads 8 \
        $PFM_PATH/debutanizer_models/basecase/vapor_liquid_simulation/simulate_vapor_liquid_submodels.jl \
        $aspen_condition_path \
        $rms_path \
        $perturb_species \
        $perturb_factor \
        $tray \
        64.0
    end=$(date +%s)
    runtime=$((end - start))
    echo "runtime: $runtime"
done
