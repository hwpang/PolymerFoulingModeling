#!/bin/bash -l

# LLsub ./submit_0.1_liquid_simulation_LLsub.sh [18,1,48] -q spot-xeon-p8
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

aspen_condition_path=$PFM_PATH/debutanizer_models/basecase/aspen_simulation/aspen_conditions.yml
rms_path=liquid_mechanism/chem_updated.rms

jobs=()
perturb_species_list=("1,3-BUTADIENE" "CYCLOPENTADIENE")
perturb_factor_list=("0.5" "0.7" "0.9" "1.0" "1.1" "1.3" "1.5" "1.7" "1.9")

for perturb_species in "${perturb_species_list[@]}"; do
    for perturb_factor in "${perturb_factor_list[@]}"; do
        jobs+=("$perturb_species $perturb_factor")
    done
done

for jobind in $(seq $LLSUB_RANK $LLSUB_SIZE ${#jobs[@]}); do
    echo "SLURM_ARRAY_TASK_ID: $LLSUB_RANK"
    echo "SLURM_ARRAY_TASK_COUNT: $LLSUB_SIZE"
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
        64.0
    end=$(date +%s)
    runtime=$((end - start))
    echo "runtime: $runtime"
done

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

chemkin_path=liquid_mechanism/chem_annotated_updated.inp
species_dict_path=liquid_mechanism/species_dictionary_updated.txt
rms_path=liquid_mechanism/chem_updated.rms
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

for jobind in $(seq $LLSUB_RANK $LLSUB_SIZE ${#jobs[@]}); do
    echo "SLURM_ARRAY_TASK_ID: $LLSUB_RANK"
    echo "SLURM_ARRAY_TASK_COUNT: $LLSUB_SIZE"
    echo "jobind: $jobind"
    echo "${jobs[$jobind]}"
    job=(${jobs[$jobind]})
    perturb_species=${job[0]}
    perturb_factor=${job[1]}
    start=$(date +%s)
    liquid_simulation_results_directory="./simulation_results/${perturb_species}_${perturb_factor}_3600.0_64.0"
    python-jl $PFM_PATH/shared/ASF_distribution/get_alpha_for_ASF_distribution.py \
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

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

conda activate rmg_py3_20230404

which python

PFM_PATH="/home/gridsan/hwpang/Software/PolymerFoulingModeling"
simulation_results_directory="simulation_results/1,3-BUTADIENE_1.0_3600.0_64.0"
aspen_condition_path=$PFM_PATH/debutanizer_models/basecase/aspen_simulation/aspen_conditions.yml
model_name="basecase_debutanizer_model"

start=$(date +%s)
python $PFM_PATH/shared/analysis/analyze_debutanizer_vapor_liquid_simulation.py \
    --model_name $model_name \
    --aspen_condition_path $aspen_condition_path \
    --simulation_results_directory $simulation_results_directory
end=$(date +%s)
runtime=$((end - start))
echo "runtime: $runtime"

num1="$((num1 + num2))"

if [[ "$LLSUB_RANK" == "$(($LLSUB_SIZE - 1))" ]]; then
    echo "============================================================"
    echo "submitting film simulation jobs"
    LLsub ./submit_0.2_film_simulation_LLsub.sh [30,24,2] -q spot-xeon-p8
    echo "============================================================"
fi
