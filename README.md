# PolymerFoulingModeling
A collection of Python and Julia scripts used for generating and simulating multi-phase chemical kinetic models for polymer fouling

## Repository structure
📦debutanizer_models
 ┣ 📂basecase # files for the base case debutanizer model
 ┃ ┣ 📂aspen_simulation # files related to process simulation
 ┃ ┃ ┣ 📜Distillation_FEDP_XN_V2_060720_v2.xlsx # output file from Aspen Plus
 ┃ ┃ ┗ 📜aspen_conditions.yml # extracted process conditions (T, P, Cmajor, etc.)
 ┃ ┣ 📂film_mechanism # mechanism files for film growth
 ┃ ┃ ┣ 📜chem_annotated_film.inp
 ┃ ┃ ┣ 📜chem_film_phase.rms
 ┃ ┃ ┣ 📜liquid_species_mapping.yml
 ┃ ┃ ┗ 📜species_dictionary_film.txt
 ┃ ┣ 📂liquid_mechanism # mechanism files for liquid phase oligomerization
 ┃ ┃ ┣ 📜chem.rms
 ┃ ┃ ┣ 📜chem_annotated.inp
 ┃ ┃ ┗ 📜species_dictionary.txt
 ┃ ┣ 📂liquid_mechanism_generation # files to generate liquid phase mechanism using RMG-Py
 ┃ ┃ ┣ 📜input.py
 ┃ ┃ ┣ 📜run_generate_input.sh
 ┃ ┃ ┗ 📜submit_input.sh
 ┃ ┗ 📂submit_scripts
 ┃ ┃ ┣ 📜submit_0.1_liquid_simulation_LLsub.sh
 ┃ ┃ ┣ 📜submit_0.2_film_simulation_LLsub.sh
 ┃ ┃ ┗ 📜submit_0_all_LLsub.sh
 ┗ 📂trace_oxygen_perturbed # files for the trace-oxygen perturbed debutanizer model
 ┃ ┣ 📂aspen_simulation
 ┃ ┃ ┣ 📜Distillation_FEDP_XN_V2_060720_v2_oxygen.xlsx
 ┃ ┃ ┗ 📜aspen_conditions.yml
 ┃ ┣ 📂film_mechanism
 ┃ ┃ ┣ 📜chem_annotated_film.inp
 ┃ ┃ ┣ 📜chem_film_phase.rms
 ┃ ┃ ┣ 📜chem_liquid_film_phase.rms
 ┃ ┃ ┣ 📜liquid_species_mapping.yml
 ┃ ┃ ┗ 📜species_dictionary_film.txt
 ┃ ┣ 📂liquid_mechanism
 ┃ ┃ ┣ 📜chem.rms
 ┃ ┃ ┣ 📜chem_annotated.inp
 ┃ ┃ ┗ 📜species_dictionary.txt
 ┃ ┣ 📂liquid_mechanism_generation
 ┃ ┃ ┗ 📜input.py
 ┃ ┗ 📂submit_scripts
 ┃ ┃ ┣ 📜submit_0.1_liquid_simulation_LLsub.sh
 ┃ ┃ ┣ 📜submit_0.2_film_simulation_LLsub.sh
 ┃ ┃ ┣ 📜submit_0.3_simulate_film_submodel_1D_LLsub.sh
 ┃ ┃ ┣ 📜submit_0_all_LLsub.sh
 ┃ ┃ ┣ 📜submit_10_generate_film_mechanism_LLsub.sh
 ┃ ┃ ┣ 📜submit_11_simulate_film_submodel_LLsub.sh
 ┃ ┃ ┣ 📜submit_12_analyze_film_simulation_and_sensitivity_LLsub.sh
 ┃ ┃ ┣ 📜submit_13_simulate_film_submodel_1D_discretization_studies_LLsub.sh
 ┃ ┃ ┣ 📜submit_14_analyze_debutanizer_film_simulation_1D_discretization_studies_LLsub.sh
 ┃ ┃ ┣ 📜submit_15_analyze_debutanizer_film_simulation_and_sensitivity_1D_LLsub.sh
 ┃ ┃ ┣ 📜submit_16_plotter.sh
 ┃ ┃ ┣ 📜submit_1_extract_and_plot_aspen_conditions_LLsub.sh
 ┃ ┃ ┣ 📜submit_2_generate_rmg_input_LLsub.sh
 ┃ ┃ ┣ 📜submit_3_generate_liquid_mechanism_LLsub.sh
 ┃ ┃ ┣ 📜submit_4_update_parameters_LLsub.sh
 ┃ ┃ ┣ 📜submit_5_time_step_studies.sh
 ┃ ┃ ┣ 📜submit_6_analyze_time_step_studies.sh
 ┃ ┃ ┣ 📜submit_7_simulate_vapor_liquid_submodels_LLsub.sh
 ┃ ┃ ┣ 📜submit_8_get_alpha_for_ASF_distribution_LLsub.sh
 ┃ ┃ ┗ 📜submit_9_analyze_vapor_liquid_simulation_LLsub.sh
 📦QCMD_cell_model
 ┣ 📂film_mechanism
 ┃ ┣ 📜chem_annotated_film.inp
 ┃ ┣ 📜chem_film_phase.rms
 ┃ ┣ 📜liquid_species_mapping.yml
 ┃ ┗ 📜species_dictionary_film.txt
 ┣ 📂liquid_mechanism
 ┃ ┣ 📜chem.rms
 ┃ ┣ 📜chem_annotated.inp
 ┃ ┗ 📜species_dictionary.txt
 ┣ 📂liquid_mechanism_generation
 ┃ ┗ 📜input.py
 ┣ 📂liquid_simulation
 ┃ ┗ 📜simulate_liquid_submodel.jl
 ┗ 📂submit_scripts
 ┃ ┣ 📜submit_0.1_simulation_LLsub.sh
 ┃ ┗ 📜submit_0_all_LLsub.sh

# Installation
- git clone RMG-Py (https://github.com/hwpang/RMG-Py), RMG-database (https://github.com/hwpang/RMG-database), and RMS (https://github.com/hwpang/ReactionMechanismSimulator.jl)
- Check out to branch `polymer_fouling_1` for RMG-Py, RMG-database, and RMS
- Install RMG-Py following http://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/installation/anacondaDeveloper.html
    - During installtion, use the following to point to your local RMS instead of the main branch on GitHub
        ```
        julia -e 'using Pkg; Pkg.develop(Pkg.PackageSpec(name="ReactionMechanismSimulator",path="/path/to/your/ReactionMechanismSimulator.jl")); Pkg.build("ReactionMechanismSimulator"); using ReactionMechanismSimulator;'
        ```
- Add these Julia packages
    ```
    julia -e 'using Pkg; Pkg.add("CSV"); Pkg.add("DataFrames"); Pkg.add("XLSX");'
    ```