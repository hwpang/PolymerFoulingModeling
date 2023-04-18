conda activate rmg_py3_20230404

which julia

PFM_PATH=/home/gridsan/hwpang/Software/PolymerFoulingModeling
rms_path=$PFM_PATH/debutanizer_models/trace_oxygen_perturbed/liquid_mechanism/chem.rms
model_name=trace_oxygen_perturbed_debutanizer_model
julia $PFM_PATH/debutanizer_models/basecase/analysis/analyze_film_simulation.jl \
        $rms_path \
        $model_name