import argparse
import yaml
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--aspen_condition_path', type=str, required=True, help='Path to the yaml file containing Aspen conditions')

    args = parser.parse_args()
    aspen_condition_path = args.aspen_condition_path
    return aspen_condition_path

aspen_condition_path = parse_arguments()

with open(aspen_condition_path, 'r') as f:
    aspen_conditions = yaml.load(f, Loader=yaml.FullLoader)

# Generate the RMG-Py input file for liquid phase mechanism generation

input_string = """
# Data sources
database(
    thermoLibraries = ['thermo_combined','multi_trays_combined','Conjugated_diene','thermo_DFT_CCSDTF12_BAC','primaryThermoLibrary', 'primaryNS', 'NitrogenCurran', 'NOx2018', 'FFCM1(-)', 'SulfurLibrary', 'SulfurGlarborgH2S'],
    reactionLibraries = ['kinetics_combined','Xu_cyclopentadiene','Conjugated_diene','Klippenstein_Glarborg2016'],
    seedMechanisms = [],
    kineticsDepositories = ['training'], 
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

# List of species

species(
    label='N-BUTANE',
    reactive=True,
    structure=SMILES("CCCC")
)

species(
    label='2-BUTENE',
    reactive=True,
    structure=SMILES("CC=CC"),
)

species(
    label='1,3-BUTADIENE',
    reactive=True,
    structure=SMILES("C=CC=C"),
)

species(
    label='CYCLOPENTADIENE',
    reactive=True,
    structure=SMILES("C1=CCC=C1"),
)

species(
    label='BENZENE',
    reactive=True,
    structure=SMILES("C1=CC=CC=C1"),
)

species(
    label='1,3-CYCLOHEXADIENE',
    reactive=True,
    structure=SMILES("C1=CCCC=C1"),
)

species(
    label='TOLUENE',
    reactive=True,
    structure=SMILES("CC1=CC=CC=C1"),
)

species(
    label='STYRENE',
    reactive=True,
    structure=SMILES("C=CC1=CC=CC=C1"),
)
"""

reactor_string = """
# tray {tray_number}
constantTVLiquidReactor(
    temperature=({temperature},'K'),
    initialConcentrations={{
{initial_concentrations}
    }},
    terminationTime=(8000,'hr'),
    residenceTime=({residence_time},'s'),
    vaporPressure = ({vapor_pressure},'Pa'),
    vaporMoleFractions={{
{vapor_mole_fractions}
    }},
)
"""

trays = range(40)
d = 2.5
h = 0.3
A = (d / 2) ** 2 * np.pi
Vliq = A * h

major_species = ['N-BUTANE', '2-BUTENE', '1,3-BUTADIENE', 'CYCLOPENTADIENE', 'BENZENE', '1,3-CYCLOHEXADIENE', 'TOLUENE', 'STYRENE']

for tray in trays:
    initial_concentrations = ''
    for species in major_species:
        concentrations = aspen_conditions['liquid_concentration'][species]
        initial_concentrations += f"        '{species}': ({concentrations[tray]}, 'mol/m^3'),"
        initial_concentrations += '\n'
    if tray+1 == len(trays):
        residence_time = Vliq/aspen_conditions["liquid_outlet_volumetric_flowrate"][tray-1] # use the inlet flow to calculate liquid residence time for reboiler instead, as most of the liquid evaporates
    else:
        residence_time = Vliq/aspen_conditions["liquid_outlet_volumetric_flowrate"][tray]
    print(initial_concentrations)
    vapor_mole_fractions = ''
    total_vapor_concentration = 0
    for species in major_species:
        concentrations = aspen_conditions['vapor_concentration'][species]
        total_vapor_concentration += concentrations[tray]
    for species in major_species:
        concentrations = aspen_conditions['vapor_concentration'][species]
        vapor_mole_fractions += f"        '{species}': {concentrations[tray]/total_vapor_concentration},"
        vapor_mole_fractions += '\n'
    print(vapor_mole_fractions)
    input_string += reactor_string.format(
        tray_number=tray+1,
        temperature=aspen_conditions["T"][tray],
        initial_concentrations=initial_concentrations,
        residence_time=residence_time,
        vapor_pressure=aspen_conditions["P"][tray],
        vapor_mole_fractions=vapor_mole_fractions,
    )

input_string += """

liquidVolumetricMassTransferCoefficientPowerLaw(
    prefactor=(14.5, "1/s"),
    diffusionCoefficientPower=1/2,
    solventViscosityPower=-1/6,
    solventDensityPower=-1/6,
)

solvation(
    solvent='benzene'
)

simulator(
    atol=1e-20,
    rtol=1e-6,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.01,
    toleranceInterruptSimulation=0.01,
    maximumEdgeSpecies=100000,
    filterReactions=True,
    toleranceBranchReactionToCore=0.01,
    branchingIndex=0.5,
    branchingRatioMax=1.0,
    maxNumObjsPerIter=10,
    toleranceReactionToCoreDeadendRadical=0.01,
)

options(
    units='si',
)

generatedSpeciesConstraints(
    #allows exceptions to the following restrictions
    allowed=['input species','seed mechanisms','reaction libraries'],
    # Constraints on generated species
    maximumRadicalElectrons = 1,
    maximumCarbonAtoms=16,
)
"""

with open('input.py', 'w') as f:
    f.write(input_string)