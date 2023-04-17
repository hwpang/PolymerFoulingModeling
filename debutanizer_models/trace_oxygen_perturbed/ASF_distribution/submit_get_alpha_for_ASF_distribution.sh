#!/bin/bash -l
#SBATCH -J get_ASF_distribution
#SBATCH -o slurm-get_ASF_distribution-%a.out
#SBATCH -n 12
#SBATCH --array=0-7

conda activate rmg_py3_20230404

which julia
which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling

chemkin_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/liquid_mechanism/chem_annotated.inp
species_dict_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/liquid_mechanism/species_dictionary.txt
rms_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/liquid_mechanism/chem.rms
aspen_condition_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/aspen_simulation/aspen_conditions_oxygen.yml
model_name="trace_oxygen_perturbed_debutanizer_model"

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
    start=`date +%s`
    liquid_simulation_results_directory="./simulation_results/${perturb_species}_${perturb_factor}_3600.0_64.0"
    python-jl $PFM_PATH/debutanizer_models/basecase/ASF_distribution/get_alpha_for_ASF_distribution.py \
                --chemkin_path $chemkin_path \
                --species_dict_path $species_dict_path \
                --rms_path $rms_path \
                --aspen_condition_path $aspen_condition_path \
                --results_directory $liquid_simulation_results_directory \
                --n_jobs 12 \
                --model_name $model_name
    end=`date +%s`
    runtime=$((end-start))
    echo "runtime: $runtime"
done
