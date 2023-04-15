#!/bin/bash -l
#SBATCH -J time_step_studies
#SBATCH -o slurm-time_step_studies-%A_%a.out
#SBATCH -n 8
#SBATCH --array=32,64,128,256,512,1024

conda activate rmg_py3_20230404

which julia
which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling

aspen_condition_path=$PFM_PATH/debutanizer_models/basecase/aspen_simulation/aspen_conditions.yml
liquid_mechanism_path=$PFM_PATH/debutanizer_models/basecase/liquid_mechanism/chem.rms

start=`date +%s`
julia --threads 8 \
        $PFM_PATH/debutanizer_models/basecase/vapor_liquid_simulation/simulate_vapor_liquid_submodels.jl \
        $aspen_condition_path \
        $liquid_mechanism_path \
        "1,3-BUTADIENE" "1.0" \
        $SLURM_ARRAY_TASK_ID
end=`date +%s`
runtime=$((end-start))
echo "runtime: $runtime"
