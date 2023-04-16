#!/bin/bash -l

conda activate rmg_py3_20230404
which julia
which python
which python-jl

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling

aspen_results_path=$PFM_PATH/debutanizer_models/DEHA_mitigated/aspen_simulation/Distillation_FEDP_XN_V2_060720_v2_oxygen_DEHA.xlsx
model_name=DEHA_mitigated_debutanizer_model

julia $PFM_PATH/debutanizer_models/basecase/aspen_simulation/extract_aspen_conditions.jl \
        $aspen_results_path \
        $model_name