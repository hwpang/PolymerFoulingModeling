#!/bin/bash -l
#SBATCH -J simulate_vapor_liquid
#SBATCH -o slurm-simulate_vapor_liquid-%a.out
#SBATCH -n 8
#SBATCH --array=0-17

conda activate rmg_py3_20230404

which julia

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling

aspen_condition_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/aspen_simulation/aspen_conditions_oxygen.yml
rms_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/liquid_mechanism/chem.rms
model_name=trace_oxygen_perturbed_debutanizer_model

jobs=()
perturb_species_list=("OXYGEN")
perturb_factor_list=("0.0" "1e-6" "1e-5" "1e-4" "1e-3" "1e-2" "1e-1" "1e0")

for perturb_species in "${perturb_species_list[@]}"
do
    for perturb_factor in "${perturb_factor_list[@]}"
    do
        jobs+=("$perturb_species $perturb_factor")
    done
done

for jobind in `seq $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT ${#jobs[@]}`
do
    echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
    echo "SLURM_ARRAY_TASK_COUNT: $SLURM_ARRAY_TASK_COUNT"
    echo "jobind: $jobind"
    echo "${jobs[$jobind]}"
    job=(${jobs[$jobind]})
    perturb_species=${job[0]}
    perturb_factor=${job[1]}
    tray=${job[2]}
    start=`date +%s`
    julia --threads 8 \
            $PFM_PATH/debutanizer_models/basecase/vapor_liquid_simulation/simulate_vapor_liquid_submodels.jl \
            $aspen_condition_path \
            $rms_path \
            $perturb_species \
            $perturb_factor \
            $tray \
            64.0 \
            $model_name
    end=`date +%s`
    runtime=$((end-start))
    echo "runtime: $runtime"
done
