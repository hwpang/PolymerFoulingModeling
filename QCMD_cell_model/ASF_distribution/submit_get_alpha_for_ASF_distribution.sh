#!/bin/bash -l
#SBATCH -J get_alpha_for_ASF_distribution
#SBATCH -o slurm-get_alpha_for_ASF_distribution-%a.out
#SBATCH -n 12
#SBATCH --array=0-4

conda activate rmg_py3_20230404

which julia
which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling

chemkin_path=$PFM_PATH/QCMD_cell_model/liquid_mechanism/chem_annotated.inp
species_dict_path=$PFM_PATH/QCMD_cell_model/liquid_mechanism/species_dictionary.txt
rms_path=$PFM_PATH/QCMD_cell_model/liquid_mechanism/chem.rms
model_name="QCMD_cell_model"

jobs=()
perturb_factor_list=("0.0" "1e-3" "1e-2" "1e-1" "1e0")

for perturb_factor in "${perturb_factor_list[@]}"
do
    jobs+=("$perturb_factor")
done

for jobind in `seq $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT ${#jobs[@]}`
do
    echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
    echo "SLURM_ARRAY_TASK_COUNT: $SLURM_ARRAY_TASK_COUNT"
    echo "jobind: $jobind"
    echo "${jobs[$jobind]}"
    job=(${jobs[$jobind]})
    perturb_factor=${job[0]}
    start=`date +%s`
    liquid_simulation_results_directory="simulation_results/O2_${perturb_factor}"
    python-jl $PFM_PATH/debutanizer_models/basecase/ASF_distribution/get_alpha_for_ASF_distribution.py \
                --chemkin_path $chemkin_path \
                --species_dict_path $species_dict_path \
                --rms_path $rms_path \
                --results_directory $liquid_simulation_results_directory \
                --n_jobs 12 \
                --model_name $model_name
    end=`date +%s`
    runtime=$((end-start))
    echo "runtime: $runtime"
done
