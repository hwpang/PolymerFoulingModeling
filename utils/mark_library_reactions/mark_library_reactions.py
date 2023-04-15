#!/usr/bin/env python
# coding: utf-8

"""
Mark the reaction family for reactions in a given liquid phase mechanism.
"""

# +
import sys

sys.path.insert(0, "/home/gridsan/hwpang/Software/RMG-Py/")

import os
import yaml
import argparse
from joblib import Parallel, delayed

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.chemkin import load_chemkin_file, save_chemkin_file, save_species_dictionary


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--chemkin_path",
        type=str,
        required=True,
        help="The path to liquid phase mechanism chemkin file",
    )
    parser.add_argument(
        "--species_dict_path",
        type=str,
        required=True,
        help="The path to species dict file",
    )
    parser.add_argument(
        "--rms_path", type=str, required=True, help="The path to rms file",
    )
    parser.add_argument(
        "--model_name",
        type=str,
        required=True,
        help="The name of the model, i.e., basecase_debutanizer_model or QCMD_cell_model.",
    )
    parser.add_argument(
        "-n",
        "--n_jobs",
        type=int,
        default=1,
        help="The number of jobs to run in parallel.",
    )
    parser.add_argument(
        "--save_directory",
        type=str,
        default=".",
        help="The directory to save the generated film phase mechanism.",
    )
    parser.add_argument(
        "--debug", action="store_true", help="Whether to run in debug mode."
    )

    args = parser.parse_args()

    chemkin_path = args.chemkin_path
    species_dict_path = args.species_dict_path
    rms_path = args.rms_path
    model_name = args.model_name
    n_jobs = args.n_jobs
    save_directory = args.save_directory
    debug = args.debug

    return (
        chemkin_path,
        species_dict_path,
        rms_path,
        n_jobs,
        model_name,
        save_directory,
        debug,
    )


# +

(
    chemkin_path,
    species_path,
    rms_path,
    n_jobs,
    model_name,
    save_directory,
    debug,
) = parse_arguments()

database = RMGDatabase()
database.load(
    path=settings["database.directory"],
    thermo_libraries=["primaryThermoLibrary"],
    kinetics_families="default",
    reaction_libraries=[],
    kinetics_depositories=["training"],
)

# +

print("loading RMG model...")
liqspcs, liqrxns = load_chemkin_file(chemkin_path, species_path)

print("loading RMS model...")
with open(rms_path) as f:
    phase_dict = yaml.load(f, Loader=yaml.Loader)
rms_spcs = phase_dict["Phases"][0]["Species"]
rms_rxns = phase_dict["Reactions"]

# Check if the RMG and RMS model match each other
if len(rms_spcs) != len(liqspcs):
    raise ValueError("The RMG and RMS model have different number of species.")
if len(rms_rxns) != len(liqrxns):
    raise ValueError("The RMG and RMS model have different number of reactions.")

print("Rename species...")

if model_name == "basecase_debutanizer_model":
    initial_monomer_labels = [
        "N-BUTANE",
        "2-BUTENE",
        "1,3-BUTADIENE",
        "CYCLOPENTADIENE",
        "BENZENE",
        "1,3-CYCLOHEXADIENE",
        "TOLUENE",
        "STYRENE",
    ]
elif model_name == "QCMD_cell_model":
    initial_monomer_labels = [
        "5-methylcyclohexadiene",
        "1-methylcyclohexadiene",
        "2-methylcyclohexadiene",
        "methylenecyclohexene",
    ]

for spc in liqspcs:
    if spc.label not in initial_monomer_labels:
        spc.label = spc.smiles

old_name_to_new_name = {}
for spc_dict in rms_spcs:
    if spc_dict["name"] not in initial_monomer_labels:
        old_name_to_new_name[spc_dict["name"]] = spc_dict["smiles"]
        spc_dict["name"] = spc_dict["smiles"]
    else:
        old_name_to_new_name[spc_dict["name"]] = spc_dict["name"]

for rxn_dict in rms_rxns:
    rxn_dict["reactants"] = [old_name_to_new_name[name] for name in rxn_dict["reactants"]]
    rxn_dict["products"] = [old_name_to_new_name[name] for name in rxn_dict["products"]]
    kinetics_dict = rxn_dict["kinetics"]
    if "efficiencies" in kinetics_dict:
        kinetics_dict["efficiencies"] = {old_name_to_new_name[name]: value for name, value in kinetics_dict["efficiencies"].items()}

# +
def get_reaction_family(rxn):
    rxns = database.kinetics.generate_reactions_from_families(
        reactants=rxn.reactants, products=rxn.products
    )
    if rxns:
        rxns[0].reactants = rxn.reactants
        rxns[0].products = rxn.products
        rxns[0].kinetics = rxn.kinetics
        rxns[0].index = rxn.index
        rxns[0].pairs = rxn.pairs
        return rxns[0]
    else:
        return rxn


needs_family_inds = [
    ind
    for ind, rxn in enumerate(liqrxns)
    if rxn.family not in database.kinetics.families
]
if debug:
    needs_family_inds = needs_family_inds[:10]

print("Marking reaction family for {} reactions...".format(len(needs_family_inds)))
marked_reactions = Parallel(n_jobs=n_jobs, verbose=5, backend="multiprocessing")(
    delayed(get_reaction_family)(liqrxns[ind]) for ind in needs_family_inds
)

for ind, marked_reaction in zip(needs_family_inds, marked_reactions):
    liqrxns[ind] = marked_reaction

# -

print("Saving liquid phase mechanism after marking reaction family...")
chemkin_save_path = os.path.join(save_directory, "chem_annotated.inp")
save_chemkin_file(chemkin_save_path, liqspcs, liqrxns)

species_save_path = os.path.join(save_directory, "species_dictionary.txt")
save_species_dictionary(species_save_path, liqspcs)

rms_save_path = os.path.join(save_directory, "chem.rms")
with open(rms_save_path, "w") as f:
    yaml.dump(phase_dict, f)

print("Loading generated liquid phase mechanism after marking reaction family...")
newliqspcs, newliqrxns = load_chemkin_file(chemkin_save_path, species_save_path)

print("Done!")
