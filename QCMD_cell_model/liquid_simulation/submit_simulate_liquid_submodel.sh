#!/bin/bash -l
#SBATCH -J simulate_liquid_submodel
#SBATCH -o slurm-simulate_liquid_submodel-%a.out
#SBATCH -n 8
#SBATCH --array=0-4

conda activate rmg_py3_20230404

which julia
which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling

rms_path=$PFM_PATH/QCMD_cell_model/liquid_mechanism/chem.rms

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
    julia --threads 8 \
            $PFM_PATH/QCMD_cell_model/liquid_simulation/simulate_liquid_submodel.jl \
            $rms_path \
            $perturb_factor
    end=`date +%s`
    runtime=$((end-start))
    echo "runtime: $runtime"
done
