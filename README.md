# PolymerFoulingModeling
A collection of Python and Julia scripts used for generating and simulating multi-phase chemical kinetic models for polymer fouling

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