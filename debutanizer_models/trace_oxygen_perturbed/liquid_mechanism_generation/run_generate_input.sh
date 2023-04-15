#!/bin/bash -l

conda activate rmg_py3_20230404
which julia
which python
which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
aspen_condition_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/aspen_simulation/aspen_conditions_oxygen.yml
python-jl $PFM_PATH/debutanizer_models/basecase/liquid_mechanism_generation/generate_input.py \
            --aspen_condition_path $aspen_condition_path \
            --model_name trace_oxygen_perturbed_debutanizer_model