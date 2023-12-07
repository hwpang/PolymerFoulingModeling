#!/bin/bash -l

# LLsub ./submit_0.3_film_simulation_1D_LLsub.sh [27,12,4] -q spot-xeon-p8
# LLsub ./submit_0.3_film_simulation_1D_LLsub.sh [8,12,4]
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

aspen_condition_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/aspen_simulation/aspen_conditions.yml
# film_mech_directory=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/film_mechanism
film_mech_directory=film_mechanism
model_name="trace_oxygen_perturbed_debutanizer_model"
time_step=32.0
final_time=3616.0

jobs=()
perturb_species_list=("OXYGEN")
perturb_factor_list=("1e0" "1e-1" "1e-2" "1e-3" "1e-4" "1e-5" "1e-6" "0.0")
# num_cell="1"
# method="same_initial_size"
num_cell="3"
method="callback"

for perturb_species in "${perturb_species_list[@]}"; do
    for perturb_factor in "${perturb_factor_list[@]}"; do
        for tray in {1..40}; do
            jobs+=("$perturb_species $perturb_factor $tray $num_cell $method")
        done
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
    num_cell=${job[3]}
    method=${job[4]}
    start=$(date +%s)
    liquid_simulation_results_path="simulation_results/${perturb_species}_${perturb_factor}_3600.0_${time_step}/simulation_vapor_liquid_yliqn_${final_time}.csv"
    julia $PFM_PATH/shared/film_simulation_1D/simulate_film_submodel_1D.jl \
        $film_mech_directory \
        $model_name \
        $liquid_simulation_results_path \
        $aspen_condition_path \
        $tray \
        $num_cell \
        $perturb_factor \
        $method \
        "false"
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

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
model_name="trace_oxygen_perturbed_debutanizer_model"
simulation_directory="simulation_results/OXYGEN_1e0_3600.0_32.0"

start=$(date +%s)
python $PFM_PATH/shared/analysis/analyze_debutanizer_film_simulation_1D.py \
    --model_name $model_name \
    --simulation_directory $simulation_directory
end=$(date +%s)
runtime=$((end - start))
echo "runtime: $runtime"

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

conda activate rmg_py3_20230404

which python

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
model_name=trace_oxygen_perturbed_debutanizer_model
all_simulation_directory="simulation_results"

start=$(date +%s)
python $PFM_PATH/shared/analysis/analyze_debutanizer_sensitivity_1D.py \
    --model_name $model_name \
    --all_simulation_directory $all_simulation_directory
end=$(date +%s)
runtime=$((end - start))
echo "runtime: $runtime"