import os
import sys
import argparse
from joblib import Parallel, delayed

from rmgpy import settings
from rmgpy.rmg.main import RMG
from rmgpy.yml import write_yml
from rmgpy.kinetics import KineticsData
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.depository import DepositoryReaction
from rmgpy.kinetics.diffusionLimited import diffusion_limiter
from rmgpy.chemkin import load_chemkin_file, save_chemkin_file, save_species_dictionary
from rmgpy.data.vaporLiquidMassTransfer import vapor_liquid_mass_transfer, liquidVolumetricMassTransferCoefficientPowerLaw

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--chemkin_path", type=str, required=True, help="Path to chemkin file",
    )
    parser.add_argument(
        "--species_dict_path", type=str, required=True, help="Path to species dictionary file",
    )
    parser.add_argument(
        "--n_jobs", type=int, default=48, help="Number of jobs to run in parallel",
    )

    args = parser.parse_args()
    chemkin_path = args.chemkin_path
    species_dict_path = args.species_dict_path
    n_jobs = args.n_jobs

    return (
        chemkin_path,
        species_dict_path,
        n_jobs,
    )

(
    chemkin_path,
    species_dict_path,
    n_jobs,
) = parse_arguments()

debug = False
include_kLA_kH = False

rmg = RMG()

rmg.database_directory = settings["database.directory"]
rmg.thermo_libraries = [
    "hwpang_fouling",
    "Conjugated_diene",
    "Klippenstein_Glarborg2016",
    "BurkeH2O2",
    "thermo_DFT_CCSDTF12_BAC",
    "DFT_QCI_thermo",
    "primaryThermoLibrary",
    "primaryNS",
    "NitrogenCurran",
    "NOx2018",
    "FFCM1(-)",
]
rmg.kinetics_families = "default"
rmg.kinetics_depositories = ["training"]
rmg.kinetics_estimator = "rate rules"
rmg.solvent = "benzene"
rmg.reaction_libraries = [
    ("hwpang_fouling", False),
    ("Xu_cyclopentadiene", False),
    ("Conjugated_diene", False),
    ("Klippenstein_Glarborg2016", False),
]

rmg.load_database()

solvent_data = rmg.database.solvation.get_solvent_data(rmg.solvent)
diffusion_limiter.enable(solvent_data, rmg.database.solvation)

if include_kLA_kH:
    prefactor = 14.5
    diffusion_coefficient_power = 1/2
    solvent_viscosity_power = -1/6
    solvent_density_power = -1/6

    liquid_volumetric_mass_transfer_coefficient_power_law = liquidVolumetricMassTransferCoefficientPowerLaw(
        prefactor=prefactor,
        diffusion_coefficient_power=diffusion_coefficient_power,
        solvent_viscosity_power=solvent_viscosity_power,
        solvent_density_power=solvent_density_power,
    )

    vapor_liquid_mass_transfer.enable(solvent_data, rmg.database.solvation, liquid_volumetric_mass_transfer_coefficient_power_law)

spcs, rxns = load_chemkin_file(
    chemkin_path,
    species_dict_path,
)

print("Loaded {0} species and {1} reactions.".format(len(spcs), len(rxns)))

model = CoreEdgeReactionModel()
model.solvent_name = rmg.solvent

def update_thermo(spc):
    self = model
    spc.thermo = None
    self.generate_thermo(spc)
    return spc.thermo, spc.liquid_volumetric_mass_transfer_coefficient_data, spc.henry_law_constant_data

def update_reaction(rxn):
    rxns = rmg.database.kinetics.generate_reactions_from_libraries(reactants=rxn.reactants, products=rxn.products,)
    if rxns:
        # rxn.family = rxns[0].family.label.replace("Kinetics Library ", "")
        # forward = rxn
        if rxns[0].kinetics is not None:
            forward = rxns[0]
        else:
            forward = rxn
        if debug:
            print("Library reaction")
            print("rxn:", forward)
            print("is forward:", forward.is_forward)
            print("reactants:", forward.reactants)
            print("kinetics:", forward.kinetics)
            print("library:", forward.library)
            print("")
    else:
        rxns = rmg.database.kinetics.generate_reactions_from_families(reactants=rxn.reactants, products=rxn.products, only_families=rxn.family)
        if rxns:
            forward = rxns[0]

            if debug:
                print("Family reaction")
                print("rxn:", forward)
                print("is forward:", forward.is_forward)
                print("reactants:", forward.reactants)
                print("kinetics", forward.kinetics)
                print("template:", forward.template)
        else:
            forward = rxn
            if debug:
                print("Weird reaction; assuming library reaction")
                print("rxn:", forward)
                print("reactants:", forward.reactants)
                print("kinetics:", forward.kinetics)
                print("")

    return forward

def apply_kinetics_to_reaction(self, reaction):
    """
    retrieve the best kinetics for the reaction and apply it towards the forward
    or reverse direction (if reverse, flip the direaction).
    """
    from rmgpy.data.rmg import get_db
    # Find the reaction kinetics
    kinetics, source, entry, is_forward = self.generate_kinetics(reaction)
    # Flip the reaction direction if the kinetics are defined in the reverse direction
    if not is_forward:
        family = get_db('kinetics').families[reaction.family]
        reaction.reactants, reaction.products = reaction.products, reaction.reactants
        reaction.pairs = [(p, r) for r, p in reaction.pairs]
        if family.own_reverse and hasattr(reaction, 'reverse'):
            if reaction.reverse:
                reaction.template = reaction.reverse.template
                # replace degeneracy
                reaction.degeneracy = reaction.reverse.degeneracy
            # We're done with the "reverse" attribute, so delete it to save a bit of memory
            reaction.reverse = None
    reaction.kinetics = kinetics
    return is_forward

def update_kinetics(forward):
    self = model
    different_direction = False
    if forward.kinetics is None:
        is_forward = apply_kinetics_to_reaction(self, forward)
        if not is_forward:
            different_direction = True

    if isinstance(forward.kinetics, KineticsData):
        forward.kinetics = forward.kinetics.to_arrhenius()
    #  correct barrier heights of estimated kinetics
    if isinstance(
        forward, (TemplateReaction, DepositoryReaction)
    ):  # i.e. not LibraryReaction
        forward.fix_barrier_height()  # also converts ArrheniusEP to Arrhenius.

    if self.pressure_dependence and forward.is_unimolecular():
        # If this is going to be run through pressure dependence code,
        # we need to make sure the barrier is positive.
        forward.fix_barrier_height(force_positive=True)
    
    if debug:
        print("forward:", forward)
        print("is forward:", forward.is_forward)
    return forward.kinetics, different_direction

print("Making new species for {0} species...".format(len(spcs)))
for spc in spcs:
    model.make_new_species(spc, check_existing=False, generate_thermo=False)

print("Updating thermo for {0} species...".format(len(spcs)))

thermos_kLAs_kHs = Parallel(n_jobs=n_jobs, verbose=5, backend="multiprocessing")(
    delayed(update_thermo)(spc) for spc in model.new_species_list
)

for spc, thermo_kLA_kH in zip(model.new_species_list, thermos_kLAs_kHs):
    thermo, liquid_volumetric_mass_transfer_coefficient_data, henry_law_constant_data = thermo_kLA_kH
    spc.thermo = thermo
    if include_kLA_kH:
        spc.liquid_volumetric_mass_transfer_coefficient_data = liquid_volumetric_mass_transfer_coefficient_data
        spc.henry_law_constant_data = henry_law_constant_data

print("Updating reactions for {0} reactions...".format(len(rxns)))

if debug:
    new_rxns = []
    for rxn in rxns[:100]:
        new_rxn = update_reaction(rxn)
        new_rxns.append(new_rxn)
else:
    new_rxns = Parallel(n_jobs=n_jobs, verbose=5, backend="multiprocessing")(
        delayed(update_reaction)(rxn) for rxn in rxns
    )

print("Making new reactions for {0} reactions...".format(len(rxns)))

for rxn in new_rxns:
    model.make_new_reaction(rxn, check_existing=False, generate_thermo=False, generate_kinetics=False, perform_cut=False,)

print("Updating kinetics for {0} reaction...".format(len(rxns)))

new_kinetics = Parallel(n_jobs=n_jobs, verbose=5, backend="multiprocessing")(
    delayed(update_kinetics)(rxn) for rxn in model.new_reaction_list
)

for rxn, kinetics_diff_dir in zip(model.new_reaction_list, new_kinetics):
    kinetics, diff_dir = kinetics_diff_dir
    if diff_dir:
        rxn.reactants, rxn.products = rxn.products, rxn.reactants
        rxn.pairs = [(p, r) for r, p in rxn.pairs]
        family = rmg.database.kinetics.families[rxn.family]
        if family.own_reverse and hasattr(rxn, 'reverse'):
            if rxn.reverse:
                rxn.template = rxn.reverse.template
                # replace degeneracy
                rxn.degeneracy = rxn.reverse.degeneracy
            # We're done with the "reverse" attribute, so delete it to save a bit of memory
            rxn.reverse = None
    rxn.kinetics = kinetics

print("Adding species and reactions to core...")
for spc in model.new_species_list:
    model.add_species_to_core(spc)

for rxn in model.new_reaction_list:
    model.add_reaction_to_core(rxn)

print("Saving updated chemkin file...")

if debug:
    suffix = "_debug"
else:
    suffix = ""

os.makedirs("liquid_mechanism", exist_ok=True)

save_chemkin_file(
    f"liquid_mechanism/chem_annotated{suffix}.inp",
    model.core.species,
    model.core.reactions,
    verbose=True,
    check_for_duplicates=True,
)

save_species_dictionary(
    f"liquid_mechanism/species_dictionary{suffix}.txt",
    model.core.species,
)

print("Saving updated rms file...")

write_yml(
    spcs=model.core.species,
    rxns=model.core.reactions,
    solvent=rmg.solvent,
    solvent_data=diffusion_limiter.solvent_data,
    path=f"liquid_mechanism/chem{suffix}.rms",
)

print("Done!")
