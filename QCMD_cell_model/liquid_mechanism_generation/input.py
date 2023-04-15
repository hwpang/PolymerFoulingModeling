# Data sources
database(
    thermoLibraries = ['ExptModelTraceO2','thermo_combined','multi_trays_combined','Conjugated_diene','Klippenstein_Glarborg2016', 'BurkeH2O2', 'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo', 'primaryThermoLibrary', 'primaryNS', 'NitrogenCurran', 'NOx2018', 'FFCM1(-)', 'SulfurLibrary', 'SulfurGlarborgH2S'],
    reactionLibraries = ['BurkeH2O2inN2'],
    seedMechanisms = [],
    kineticsDepositories = ['training'], 
    kineticsFamilies = ['default'],
    kineticsEstimator = 'rate rules',
)

# List of species

species(
    label='5-methylcyclohexadiene',
    reactive=True,
    structure=SMILES("CC1C=CC=CC1")
)

species(
    label='1-methylcyclohexadiene',
    reactive=True,
    structure=SMILES("CC1=CC=CCC1")
)

species(
    label='2-methylcyclohexadiene',
    reactive=True,
    structure=SMILES("CC1=CCCC=C1")
)

species(
    label='methylenecyclohexene',
    reactive=True,
    structure=SMILES("C=C1C=CCCC1")
)

species(
    label='benzene',
    reactive=False,
    structure=SMILES("c1ccccc1")
)

species(
    label='O2',
    reactive=True,
    structure=SMILES("[O][O]")
)

# Reaction system 1
constantTVLiquidReactor(
    temperature=(90 + 273.15,'K'),
    initialConcentrations={
        "benzene": (7850.467289719625, 'mol/m^3'),
        "5-methylcyclohexadiene": (73.60594796, 'mol/m^3'),
        "1-methylcyclohexadiene": (1153.1598513, 'mol/m^3'),
        "2-methylcyclohexadiene": (1177.69516729, 'mol/m^3'),
        "methylenecyclohexene": (49.07063197, 'mol/m^3'),
    },
    terminationTime=(2,'hr'),
    residenceTime=(16,'s'),
)

# Reaction system 2
constantTVLiquidReactor(
    temperature=(90 + 273.15,'K'),
    initialConcentrations={
        "benzene": (7850.467289719625, 'mol/m^3'),
        "5-methylcyclohexadiene": (73.60594796, 'mol/m^3'),
        "1-methylcyclohexadiene": (1153.1598513, 'mol/m^3'),
        "2-methylcyclohexadiene": (1177.69516729, 'mol/m^3'),
        "methylenecyclohexene": (49.07063197, 'mol/m^3'),
        "O2": (100.0, 'mol/m^3'), #1e-3 M O2
    },
    terminationTime=(2,'hr'),
    residenceTime=(16,'s'),
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
    maximumRadicalElectrons = 2,
    maximumCarbonAtoms=18,
)

