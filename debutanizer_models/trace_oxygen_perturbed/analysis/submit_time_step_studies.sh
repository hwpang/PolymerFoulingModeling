#!/bin/bash -l
#SBATCH -J time_step_studies
#SBATCH -o slurm-time_step_studies-%a.out
#SBATCH -n 8
#SBATCH --array=32,64,128,256,512,1024

conda activate rmg_py3_20230404
which julia
which python

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling

aspen_condition_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/aspen_simulation/aspen_conditions_oxygen.yml
liquid_mechanism_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/liquid_mechanism/chem.rms

start=`date +%s`
julia --project=$project_path \
        --threads 8 \
        $PFM_PATH/debutanizer_models/basecase/vapor_liquid_simulation/simulate_vapor_liquid_submodels.jl \
        $aspen_condition_path \
        $liquid_mechanism_path \
        "OXYGEN" "1.0" \
        $SLURM_ARRAY_TASK_ID \
        "trace_oxygen_perturbed_debutanizer_model"
end=`date +%s`
runtime=$((end-start))
echo "runtime: $runtime"
