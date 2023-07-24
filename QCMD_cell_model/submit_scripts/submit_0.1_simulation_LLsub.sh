#!/bin/bash -l

# LLsub ./submit_0_simulation_LLsub.sh [5,1,48] -q spot-xeon-p8
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

rms_path=liquid_mechanism/chem_updated.rms

jobs=()
perturb_factor_list=("0.0" "1e-3" "1e-2" "1e-1" "1e0")

for perturb_factor in "${perturb_factor_list[@]}"; do
    jobs+=("$perturb_factor")
done

for jobind in $(seq $LLSUB_RANK $LLSUB_SIZE ${#jobs[@]}); do
    echo "LLSUB_RANK: $LLSUB_RANK"
    echo "LLSUB_SIZE: $LLSUB_SIZE"
    echo "jobind: $jobind"
    echo "${jobs[$jobind]}"
    job=(${jobs[$jobind]})
    perturb_factor=${job[0]}
    start=$(date +%s)
    julia --threads 8 \
        $PFM_PATH/QCMD_cell_model/liquid_simulation/simulate_liquid_submodel.jl \
        $rms_path \
        $perturb_factor
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
which python

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
model_name=QCMD_cell_model
all_simulation_directory="simulation_results"

start=$(date +%s)
python $PFM_PATH/shared/analysis/analyze_QCMD_liquid_simulation.py \
    --model_name $model_name \
    --all_simulation_directory $all_simulation_directory
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

which julia
which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling

chemkin_path=liquid_mechanism/chem_annotated_updated.inp
species_dict_path=liquid_mechanism/species_dictionary_updated.txt
rms_path=liquid_mechanism/chem_updated.rms
model_name="QCMD_cell_model"

jobs=()
perturb_factor_list=("0.0" "1e-3" "1e-2" "1e-1" "1e0")

for perturb_factor in "${perturb_factor_list[@]}"; do
    jobs+=("$perturb_factor")
done

for jobind in $(seq $LLSUB_RANK $LLSUB_SIZE ${#jobs[@]}); do
    echo "LLSUB_RANK: $LLSUB_RANK"
    echo "LLSUB_SIZE: $LLSUB_SIZE"
    echo "jobind: $jobind"
    echo "${jobs[$jobind]}"
    job=(${jobs[$jobind]})
    perturb_factor=${job[0]}
    start=$(date +%s)
    liquid_simulation_results_directory="simulation_results/O2_${perturb_factor}"
    python-jl $PFM_PATH/shared/ASF_distribution/get_alpha_for_ASF_distribution.py \
        --chemkin_path $chemkin_path \
        --species_dict_path $species_dict_path \
        --rms_path $rms_path \
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

which julia
which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling

rms_mech_dir=film_mechanism
model_name="QCMD_cell_model"

jobs=()
perturb_factor_list=("0.0" "1e-3" "1e-2" "1e-1" "1e0")

for perturb_factor in "${perturb_factor_list[@]}"; do
    jobs+=("$perturb_factor")
done

for jobind in $(seq $LLSUB_RANK $LLSUB_SIZE ${#jobs[@]}); do
    echo "LLSUB_RANK: $LLSUB_RANK"
    echo "LLSUB_SIZE: $LLSUB_SIZE"
    echo "jobind: $jobind"
    job=(${jobs[$jobind]})
    perturb_factor=${job[0]}
    echo "perturb_factor: $perturb_factor"
    liquid_simulation_results_path="simulation_results/O2_$perturb_factor/simulation_liquid_1.csv"
    echo "liquid_simulation_results_path: $liquid_simulation_results_path"
    start=$(date +%s)
    julia $PFM_PATH/shared/film_simulation/simulate_film_submodel.jl \
        $rms_mech_dir \
        $model_name \
        $liquid_simulation_results_path
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
which python

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
model_name=QCMD_cell_model
all_simulation_directory="simulation_results"

start=$(date +%s)
python $PFM_PATH/shared/analysis/analyze_QCMD_film_simulation.py \
    --model_name $model_name \
    --all_simulation_directory $all_simulation_directory
end=$(date +%s)
runtime=$((end - start))
echo "runtime: $runtime"
