#!/usr/bin/env python
# coding: utf-8

"""
Generate a fragment based film growth mechanism in .rms for a given liquid phase mechanism in chemkin format.
"""

import sys

sys.path.insert(0, "/home/gridsan/hwpang/Software/RMG-Py/")

import time
import os
import yaml
import random
import argparse
from joblib import Parallel, delayed

from rmgpy import settings
from rmgpy.rmg.main import RMG
from rmgpy.yml import write_yml
from rmgpy.species import Species
from rmgpy.molecule.group import Group
from rmgpy.kinetics import KineticsData
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.kinetics.diffusionLimited import diffusion_limiter
from rmgpy.data.kinetics.depository import DepositoryReaction
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
        help="The path including species dict file",
    )
    parser.add_argument(
        "--model_name",
        type=str,
        required=True,
        help="The name of the model, i.e., basecase_debutanizer_model or QCMD_cell_model.",
    )
    parser.add_argument(
        "--n_jobs",
        type=int,
        default=1,
        help="The number of jobs to run in parallel.",
    )
    parser.add_argument(
        "--debug", action="store_true", help="Whether to run in debug mode."
    )
    parser.add_argument(
        "--save_directory",
        type=str,
        default="film_mechanism",
        help="The directory to save the generated film phase mechanism.",
    )

    args = parser.parse_args()

    chemkin_path = args.chemkin_path
    species_dict_path = args.species_dict_path
    n_jobs = args.n_jobs
    model_name = args.model_name
    debug = args.debug
    save_directory = args.save_directory

    return chemkin_path, species_dict_path, n_jobs, model_name, debug, save_directory


start = time.time()

print("Parsing arguments...")
(
    chemkin_path,
    species_dict_path,
    n_jobs,
    model_name,
    debug,
    save_directory,
) = parse_arguments()

print(f"Chemkin path: {chemkin_path}")
print(f"Species dict path: {species_dict_path}")
print(f"Using {n_jobs} workers...")
print(f"Generating film phase submodel for {model_name}...")
print(f"Debug mode: {debug}")

project_path = os.path.dirname(chemkin_path)

if model_name in [
    "basecase_debutanizer_model",
    "trace_oxygen_perturbed_debutanizer_model",
]:
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
    if model_name == "trace_oxygen_perturbed_debutanizer_model":
        include_oxygen = True
    else:
        include_oxygen = False
elif model_name == "QCMD_cell_model":
    initial_monomer_labels = [
        "5-methylcyclohexadiene",
        "1-methylcyclohexadiene",
        "2-methylcyclohexadiene",
        "methylenecyclohexene",
    ]
    include_oxygen = True

print("Loading liquid phase submodel...")

liqspcs, liqrxns = load_chemkin_file(
    os.path.join(chemkin_path),
    dictionary_path=os.path.join(species_dict_path),
)

print("Initializing RMG object for generating film phase submodel...")
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

print("Labeling reaction families...")


def get_reaction_family(rxn):
    if rxn.family in rmg.database.kinetics.families:
        return rxn.family
    else:
        rxns = rmg.database.kinetics.generate_reactions_from_families(
            reactants=rxn.reactants, products=rxn.products
        )
        if rxns:
            print(rxns[0].family)
            return rxns[0].family
        else:
            return rxn.family


families = Parallel(n_jobs=n_jobs, verbose=5, backend="multiprocessing")(
    delayed(get_reaction_family)(rxn) for rxn in liqrxns
)

for rxn, family in zip(liqrxns, families):
    rxn.family = family

# create helper dictionaries
liq_smiles_to_spc_dict = dict()
liq_label_to_spc_dict = dict()
liq_smiles_to_label = dict()
for spc in liqspcs:
    for mol in spc.molecule:
        liq_smiles_to_spc_dict[mol.smiles] = spc
        liq_smiles_to_label[mol.smiles] = spc.label
    liq_label_to_spc_dict[spc.label] = spc

if debug:
    liqspcs = liqspcs[:50]
    liqrxns = liqrxns[:100]

# Useful structure definitions

group_dict = dict()
group_dict["allylic_CH"] = Group().from_adjacency_list(
    """
1  C u0 p0 c0 {2,D}
2  C u0 p0 c0 {1,D} {3,S}
3 * Cs u0 p0 c0 {2,S} {4,S}
4  H u0 p0 c0 {3,S}
"""
)
group_dict["allylic_C."] = Group().from_adjacency_list(
    """
1  C u0 p0 c0 {2,D}
2  C u0 p0 c0 {1,D} {3,S}
3 * Cs u1 p0 c0 {2,S}
"""
)
group_dict["alkyl_C."] = Group().from_adjacency_list(
    """
1 Cs u1 p0 c0
"""
)
group_dict["phenyl_C."] = Group().from_adjacency_list(
    """
1 Cs u1 p0 c0 {2,S}
2 Cb u0 p0 c0 {1,S}
"""
)
group_dict["conjugated_diene"] = Group().from_adjacency_list(
    """
1  C u0 p0 c0 {2,D}
2  C u0 p0 c0 {1,D} {3,S}
3  C u0 p0 c0 {2,S} {4,D}
4  C u0 p0 c0 {3,D}
"""
)
group_dict["carbon_double_bond"] = Group().from_adjacency_list(
    """
1  C u0 p0 c0 {2,D}
2  C u0 p0 c0 {1,D}
"""
)
if include_oxygen:
    group_dict["COO."] = Group().from_adjacency_list(
        """
1  C u0 p0 c0 {2,S}
2  O u0 p2 c0 {1,S} {3,S}
3  O u1 p2 c0 {2,S}
"""
    )
    group_dict["COOC"] = Group().from_adjacency_list(
        """
1  C u0 p0 c0 {2,S}
2  O u0 p2 c0 {1,S} {3,S}
3  O u0 p2 c0 {2,S} {4,S}
4  C u0 p0 c0 {3,S}
"""
    )
    group_dict["COOH"] = Group().from_adjacency_list(
        """
1  C u0 p0 c0 {2,S}
2  O u0 p2 c0 {1,S} {3,S}
3  O u0 p2 c0 {2,S} {4,S}
4  H u0 p0 c0 {3,S}
"""
    )
    group_dict["CO."] = Group().from_adjacency_list(
        """
1  C u0 p0 c0 {2,S}
2  O u1 p2 c0 {1,S}
"""
    )

spc = Species().from_smiles("CC=C")
assert spc.molecule[0].is_subgraph_isomorphic(group_dict["allylic_CH"])
spc = Species().from_smiles("CC=C")
assert spc.molecule[0].is_subgraph_isomorphic(group_dict["allylic_CH"])
spc = Species().from_smiles("C=CC=C")
assert spc.molecule[0].is_subgraph_isomorphic(group_dict["allylic_CH"]) == False
spc = Species().from_smiles("C=C[CH2]")
assert spc.molecule[0].is_subgraph_isomorphic(group_dict["allylic_C."])
if include_oxygen:
    spc = Species().from_smiles("CO[O]")
    assert spc.molecule[0].is_subgraph_isomorphic(group_dict["COO."])


def calculate_subgraph_isomorphisms(mol, group):
    if group == "allylic_CH":
        sssr = mol.get_smallest_set_of_smallest_rings()
        subgraph_isomorphisms = mol.find_subgraph_isomorphisms(group_dict["allylic_CH"])
        count = len(subgraph_isomorphisms)
        for subgraph_isomorphism in subgraph_isomorphisms:
            for atom, group_atom in subgraph_isomorphism.items():
                if group_atom.label == "*":
                    break
            if sum(atom in ring for ring in sssr) > 1:
                count -= 1
        return count
    elif group == "allylic_C.":
        sssr = mol.get_smallest_set_of_smallest_rings()
        subgraph_isomorphisms = mol.find_subgraph_isomorphisms(group_dict["allylic_C."])
        count = len(subgraph_isomorphisms)
        for subgraph_isomorphism in subgraph_isomorphisms:
            for atom, group_atom in subgraph_isomorphism.items():
                if group_atom.label == "*":
                    break
            if sum(atom in ring for ring in sssr) > 1:
                count -= 1
        return count > 0
    elif group == "alkyl_C.":
        return (
            mol.is_subgraph_isomorphic(group_dict["alkyl_C."])
            and not mol.is_subgraph_isomorphic(group_dict["phenyl_C."])
            and not mol.is_subgraph_isomorphic(group_dict["allylic_C."])
        )
    elif group == "conjugated_diene":
        return len(mol.find_subgraph_isomorphisms(group_dict["conjugated_diene"])) // 2
    elif group == "carbon_double_bond":
        return (
            len(mol.find_subgraph_isomorphisms(group_dict["carbon_double_bond"])) // 2
        )
    elif group == "COO.":
        return len(mol.find_subgraph_isomorphisms(group_dict["COO."]))
    elif group == "COOC":
        return len(mol.find_subgraph_isomorphisms(group_dict["COOC"])) // 2
    elif group == "COOH":
        return len(mol.find_subgraph_isomorphisms(group_dict["COOH"]))
    elif group == "CO.":
        return len(mol.find_subgraph_isomorphisms(group_dict["CO."]))


# Dummy solid phase species place holder
fragment_species = dict()
fragment_species["CD"] = Species().from_smiles("CC=CC=CC")
fragment_species["CDB"] = Species().from_smiles("CCC=CCC")
fragment_species["AH"] = Species().from_smiles("C=C(C)C")
fragment_species["AR"] = Species().from_smiles("C=C(C)[CH2]")
fragment_species["KR"] = Species().from_smiles("CC(C)[CH2]")
fragment_species["inert(S)"] = Species().from_smiles("C" * 21)
if include_oxygen:
    fragment_species["CP"] = Species().from_smiles("CC(C)COOCC(C)C")
    fragment_species["HP"] = Species().from_smiles("CC(C)COO")
    fragment_species["PR"] = Species().from_smiles("CC(C)CO[O]")
    fragment_species["OR"] = Species().from_smiles("CC(C)C[O]")
    fragment_species["OH"] = Species().from_smiles("CC(C)CO")
for spc in fragment_species.values():
    spc.generate_resonance_structures()

fragment_species_implicit_structures = dict()
fragment_species_implicit_structures["CD"] = ["AH"] * calculate_subgraph_isomorphisms(
    fragment_species["CD"].molecule[0], "allylic_CH"
)
fragment_species_implicit_structures["CDB"] = ["AH"] * calculate_subgraph_isomorphisms(
    fragment_species["CDB"].molecule[0], "allylic_CH"
)
fragment_species_implicit_structures["AH"] = ["CDB"] + ["AH"] * (
    calculate_subgraph_isomorphisms(fragment_species["CDB"].molecule[0], "allylic_CH")
    - 1
)
fragment_species_implicit_structures["AR"] = ["CDB"] + [
    "AH"
] * calculate_subgraph_isomorphisms(fragment_species["CDB"].molecule[0], "allylic_CH")

fragment_radical_labels = ["AR", "KR"]
if include_oxygen:
    fragment_radical_labels += ["PR", "OR"]
fragment_radical_smiles = [
    mol.smiles
    for label in fragment_radical_labels
    for mol in fragment_species[label].molecule
]

fragment_H_labels = ["AH"]
if include_oxygen:
    fragment_H_labels += ["HP"]
fragment_H_smiles = [
    mol.smiles
    for label in fragment_H_labels
    for mol in fragment_species[label].molecule
]

for label, spc in fragment_species.items():
    spc.label = label

print("Generating film phase submodel...")

film_phase_model = CoreEdgeReactionModel()
film_phase_model.solvent_name = "benzene"

for spc in fragment_species.values():
    film_phase_model.make_new_species(spc, generate_thermo=False)


def generate_liq_film_reactions(spc):
    liq_film_reactions = []

    if (
        calculate_subgraph_isomorphisms(spc.molecule[0], "conjugated_diene") > 0
        and not spc.molecule[0].is_radical()
    ):
        liq_film_reactions += rmg.database.kinetics.generate_reactions_from_families(
            reactants=[fragment_species["CDB"], spc],
            only_families=["Diels_alder_addition"],
        )
        liq_film_reactions += rmg.database.kinetics.generate_reactions_from_families(
            reactants=[fragment_species["AR"], spc],
            only_families=["R_Addition_MultipleBond"],
        )
        liq_film_reactions += rmg.database.kinetics.generate_reactions_from_families(
            reactants=[fragment_species["KR"], spc],
            only_families=["R_Addition_MultipleBond"],
        )
        if include_oxygen:
            liq_film_reactions += (
                rmg.database.kinetics.generate_reactions_from_families(
                    reactants=[fragment_species["PR"], spc],
                    only_families=["R_Addition_MultipleBond"],
                )
            )
            liq_film_reactions += (
                rmg.database.kinetics.generate_reactions_from_families(
                    reactants=[fragment_species["OR"], spc],
                    only_families=["R_Addition_MultipleBond"],
                )
            )

    if (
        calculate_subgraph_isomorphisms(spc.molecule[0], "conjugated_diene") > 0
        and not spc.molecule[0].is_radical()
    ):
        liq_film_reactions += rmg.database.kinetics.generate_reactions_from_families(
            reactants=[fragment_species["CD"], spc],
            only_families=["Diels_alder_addition"],
        )

    # if calculate_subgraph_isomorphisms(spc.molecule[0], "carbon_double_bond") > 0 and not spc.molecule[0].is_radical():
    #     liq_film_reactions += rmg.database.kinetics.generate_reactions_from_families(reactants=[fragment_species["KR"],spc],only_families=["R_Addition_MultipleBond"])

    if spc.multiplicity == 2:
        liq_film_reactions += rmg.database.kinetics.generate_reactions_from_families(
            reactants=[fragment_species["CDB"], spc],
            only_families=["R_Addition_MultipleBond"],
        )

        # rxns = rmg.database.kinetics.generate_reactions_from_families(reactants=[fragment_species["AH"],spc],only_families=["H_Abstraction"])
        # for reaction in rxns:
        #     if any(mol.smiles == fragment_species["AR"].smiles for new_allylic_C in reaction.reactants+reaction.products for mol in new_allylic_C.molecule):
        #         if any(new_radical.smiles in liq_smiles_to_spc_dict and new_radical.multiplicity==1 for new_radical in reaction.reactants+reaction.products):
        #             liq_film_reactions.append(reaction)

        # if include_oxygen:
        #     rxns = rmg.database.kinetics.generate_reactions_from_families(reactants=[fragment_species["HP"],spc],only_families=["H_Abstraction"])
        #     for reaction in rxns:
        #         if any(mol.smiles == fragment_species["PR"].smiles for new_COO in reaction.reactants+reaction.products for mol in new_COO.molecule):
        #             if any(new_radical.smiles in liq_smiles_to_spc_dict and new_radical.multiplicity==1 for new_radical in reaction.reactants+reaction.products):
        #                 liq_film_reactions.append(reaction)

    if (
        calculate_subgraph_isomorphisms(spc.molecule[0], "allylic_CH") > 0
        and not spc.molecule[0].is_radical()
    ):
        for fragment_radical_label in fragment_radical_labels:
            rxns = rmg.database.kinetics.generate_reactions_from_families(
                reactants=[fragment_species[fragment_radical_label], spc],
                only_families=["H_Abstraction"],
            )
            for reaction in rxns:
                for new_allylic_C in reaction.reactants + reaction.products:
                    if (
                        new_allylic_C.multiplicity == 2
                        and calculate_subgraph_isomorphisms(
                            new_allylic_C.molecule[0], "allylic_C."
                        )
                        > 0
                    ):
                        if all(
                            mol.smiles != fragment_species["AR"].smiles
                            for mol in new_allylic_C.molecule
                        ):
                            if any(
                                mol.smiles in liq_smiles_to_spc_dict
                                for mol in new_allylic_C.molecule
                            ):
                                liq_film_reactions.append(reaction)

        liq_film_reactions += rmg.database.kinetics.generate_reactions_from_families(
            reactants=[fragment_species["CDB"], spc], only_families=["Retroene"]
        )

    if include_oxygen:
        # if calculate_subgraph_isomorphisms(spc.molecule[0], "COOH") > 0 and not spc.molecule[0].is_radical():

        #     for fragment_radical_label in fragment_radical_labels:
        #         rxns = rmg.database.kinetics.generate_reactions_from_families(reactants=[fragment_species[fragment_radical_label],spc],only_families=["H_Abstraction"])
        #         for reaction in rxns:
        #             for new_COO in reaction.reactants+reaction.products:
        #                 if new_COO.multiplicity == 2 and calculate_subgraph_isomorphisms(new_COO.molecule[0], "COO.") > 0:
        #                     if all(mol.smiles != fragment_species["PR"].smiles for mol in new_allylic_C.molecule):
        #                         if any(mol.smiles in liq_smiles_to_spc_dict for mol in new_COO.molecule):
        #                             liq_film_reactions.append(reaction)

        if spc.smiles == "[O][O]":
            liq_film_reactions += (
                rmg.database.kinetics.generate_reactions_from_families(
                    reactants=[fragment_species["AR"], spc],
                    only_families=["R_Recombination"],
                )
            )
            liq_film_reactions += (
                rmg.database.kinetics.generate_reactions_from_families(
                    reactants=[fragment_species["KR"], spc],
                    only_families=["R_Recombination"],
                )
            )

        if spc.smiles == "[O]O":
            liq_film_reactions += (
                rmg.database.kinetics.generate_reactions_from_families(
                    reactants=[fragment_species["CDB"], spc],
                    only_families=["HO2_Elimination_from_PeroxyRadical"],
                )
            )

        if spc.smiles == "[OH]":
            liq_film_reactions += (
                rmg.database.kinetics.generate_reactions_from_families(
                    reactants=[fragment_species["OR"], spc],
                    only_families=["R_Recombination"],
                )
            )
    return liq_film_reactions


start_1 = time.time()

liq_film_reactions = Parallel(n_jobs=n_jobs, verbose=5, backend="multiprocessing")(
    delayed(generate_liq_film_reactions)(spc) for spc in liqspcs
)
liq_film_reactions = [r for rxns in liq_film_reactions for r in rxns]
if include_oxygen:
    liq_film_reactions += rmg.database.kinetics.generate_reactions_from_families(
        reactants=[fragment_species["OR"], fragment_species["OR"]],
        only_families=["R_Recombination"],
    )

end_1 = time.time()
print("Time to generate film reactions: {}".format(end_1 - start_1))

# Dummy liquid phase oligomer place holder
print("Generating oligomer reactions...")
start_1 = time.time()
oligomer_species = dict()
oligomer_species["C=C(L,oligomer)"] = Species().from_smiles("CC(C)=CC")
oligomer_species["allylic_CH(L,oligomer)"] = Species().from_smiles("CC(C)=CCC")
oligomer_species["allylic_C.(L,oligomer)"] = Species().from_smiles("CC([CH2])=CCC")
oligomer_species["alkyl_C.(L,oligomer)"] = Species().from_smiles("CC([CH2])CCC")
if include_oxygen:
    oligomer_species["COO.(L,oligomer)"] = Species().from_smiles("CC(CO[O])CCC")
    oligomer_species["COOH(L,oligomer)"] = Species().from_smiles("CC(COO)CCC")
for spc in oligomer_species.values():
    spc.generate_resonance_structures()

oligomer_radical_labels = ["allylic_C.(L,oligomer)", "alkyl_C.(L,oligomer)"]
if include_oxygen:
    oligomer_radical_labels.append("COO.(L,oligomer)")
oligomer_radical_smiles = [
    mol.smiles
    for label in oligomer_radical_labels
    for mol in oligomer_species[label].molecule
]

oligomer_H_labels = ["allylic_CH(L,oligomer)"]
if include_oxygen:
    oligomer_H_labels.append("COOH(L,oligomer)")

for label, spc in oligomer_species.items():
    spc.label = label

for spc in oligomer_species.values():
    film_phase_model.make_new_species(spc, generate_thermo=False)

liq_film_reactions += rmg.database.kinetics.generate_reactions_from_families(
    reactants=[fragment_species["CDB"], oligomer_species["allylic_CH(L,oligomer)"]],
    only_families=["Retroene"],
)

for oligomer_radical_label in oligomer_radical_labels:
    liq_film_reactions += rmg.database.kinetics.generate_reactions_from_families(
        reactants=[fragment_species["CDB"], oligomer_species[oligomer_radical_label]],
        only_families=["R_Addition_MultipleBond"],
    )
    liq_film_reactions += rmg.database.kinetics.generate_reactions_from_families(
        reactants=[fragment_species["CD"], oligomer_species[oligomer_radical_label]],
        only_families=["R_Addition_MultipleBond"],
    )

liq_film_reactions += rmg.database.kinetics.generate_reactions_from_families(
    reactants=[fragment_species["KR"], oligomer_species["C=C(L,oligomer)"]],
    only_families=["R_Addition_MultipleBond"],
)

liq_film_reactions += rmg.database.kinetics.generate_reactions_from_families(
    reactants=[fragment_species["AR"], oligomer_species["C=C(L,oligomer)"]],
    only_families=["R_Addition_MultipleBond"],
)

liq_film_reactions += rmg.database.kinetics.generate_reactions_from_families(
    reactants=[fragment_species["CD"], oligomer_species["C=C(L,oligomer)"]],
    only_families=["Diels_alder_addition"],
)

for fragment_radical_label in fragment_radical_labels:
    for oligomer_H_label in oligomer_H_labels:
        rxns = rmg.database.kinetics.generate_reactions_from_families(
            reactants=[
                fragment_species[fragment_radical_label],
                oligomer_species[oligomer_H_label],
            ],
            only_families=["H_Abstraction"],
        )
        for reaction in rxns:
            if any(
                spc.smiles in oligomer_radical_smiles
                for spc in reaction.reactants + reaction.products
            ):
                if fragment_radical_label == "AR":
                    if any(
                        mol.smiles == fragment_species["AH"].smiles
                        for spc in reaction.reactants + reaction.products
                        for mol in spc.molecule
                    ):
                        liq_film_reactions.append(reaction)
                elif fragment_radical_label == "PR":
                    if any(
                        mol.smiles == fragment_species["HP"].smiles
                        for spc in reaction.reactants + reaction.products
                        for mol in spc.molecule
                    ):
                        liq_film_reactions.append(reaction)
                else:  # KR, OR
                    liq_film_reactions.append(reaction)

# for oligomer_radical_label in oligomer_radical_labels:
#     for fragment_H_label in fragment_H_labels:
#         rxns = rmg.database.kinetics.generate_reactions_from_families(reactants=[fragment_species[fragment_H_label],oligomer_species[oligomer_radical_label]],only_families=["H_Abstraction"])
#         for reaction in rxns:
#             if any(spc.smiles in fragment_radical_smiles for spc in reaction.reactants+reaction.products):
#                 if oligomer_radical_label == "allylic_C.(L,oligomer)":
#                     if any(spc.smiles == oligomer_species["allylic_CH(L,oligomer)"].smiles for spc in reaction.reactants+reaction.products):
#                         liq_film_reactions.append(reaction)
#                 elif oligomer_radical_label == "COO.(L,oligomer)":
#                     if any(spc.smiles == oligomer_species["COOH(L,oligomer)"].smiles for spc in reaction.reactants+reaction.products):
#                         liq_film_reactions.append(reaction)
#                 else: #COO.(L,oligomer), alkyl_C.(L,oligomer)
#                     liq_film_reactions.append(reaction)

end_1 = time.time()
print("Time to generate oligomer reactions: {}".format(end_1 - start_1))

print("Getting library kinetics...")


def get_library_reactions(reaction):
    rxns = rmg.database.kinetics.generate_reactions_from_libraries(
        reactants=reaction.reactants, products=reaction.products
    )
    if rxns:
        print(reaction.family)
        rxn = rxns[0]
        rxn.family = reaction.family
        return rxn
    else:
        return reaction


start_1 = time.time()
new_liq_film_reactions = Parallel(n_jobs=n_jobs, verbose=5, backend="multiprocessing")(
    delayed(get_library_reactions)(reaction) for reaction in liq_film_reactions
)
end_1 = time.time()
print("Time to get library reactions: {}".format(end_1 - start_1))

print("Making new reactions...")
start_1 = time.time()
for reaction, new_reaction in zip(liq_film_reactions, new_liq_film_reactions):
    new_reaction.family = reaction.family

    film_phase_model.make_new_reaction(
        new_reaction,
        generate_kinetics=False,
        generate_thermo=False,
        perform_cut=False,
    )
end_1 = time.time()
print("Time to make new reactions: {}".format(end_1 - start_1))

for spc in liqspcs:
    film_phase_model.make_new_species(spc, generate_thermo=False)


def generate_thermo(spc):
    film_phase_model.generate_thermo(spc)
    return spc.thermo


print("Generate thermos...")
start_1 = time.time()
thermos = Parallel(n_jobs=n_jobs, verbose=5, backend="multiprocessing")(
    delayed(generate_thermo)(spc) for spc in film_phase_model.new_species_list
)
for spc, thermo in zip(film_phase_model.new_species_list, thermos):
    spc.thermo = thermo
end_1 = time.time()
print("Time to generate thermos: {}".format(end_1 - start_1))


def generate_kinetics(forward):
    self = film_phase_model
    if forward.kinetics is None:
        self.apply_kinetics_to_reaction(forward)

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
    return forward.kinetics


print("Generate kinetics...")
start_1 = time.time()
kinetics = Parallel(n_jobs=n_jobs, verbose=5, backend="multiprocessing")(
    delayed(generate_kinetics)(forward)
    for forward in film_phase_model.new_reaction_list
)
for forward, kinetics in zip(film_phase_model.new_reaction_list, kinetics):
    forward.kinetics = kinetics
end_1 = time.time()
print("Time to generate kinetics: {}".format(end_1 - start_1))

print("Adding new species and reactions to core...")
for spc in film_phase_model.new_species_list:
    film_phase_model.add_species_to_core(spc)

for reaction in film_phase_model.new_reaction_list:
    if reaction.family == "R_Addition_MultipleBond" and len(reaction.products) != 1:
        continue
    film_phase_model.add_reaction_to_core(reaction)

# for label, spc in fragment_species.items():
#     for spec in film_phase_model.core.species:
#         if any(spc.smiles == mol.smiles for mol in spec.molecule):
#             spec.label = label
#             fragment_species[label] = spec

# for label, spc in oligomer_species.items():
#     for spec in film_phase_model.core.species:
#         if any(spc.smiles == mol.smiles for mol in spec.molecule):
#             spec.label = label
#             oligomer_species[label] = spec

for spc in film_phase_model.core.species:
    for mol in spc.molecule:
        if mol.smiles in liq_smiles_to_label:
            spc.label = liq_smiles_to_label[mol.smiles] + "(L)"


print("Saving film phase submodel...")
print(
    f"{len(film_phase_model.core.species)} species and {len(film_phase_model.core.reactions)} reactions in film phase submodel"
)

film_phase_model.mark_chemkin_duplicates()

if not os.path.exists(save_directory):
    os.makedirs(save_directory)

write_yml(
    spcs=film_phase_model.core.species,
    rxns=film_phase_model.core.reactions,
    solvent=film_phase_model.solvent_name,
    solvent_data=diffusion_limiter.solvent_data,
    path=os.path.join(save_directory, f"chem_film_phase.rms"),
)

save_chemkin_file(
    os.path.join(save_directory, "chem_annotated_film.inp"),
    film_phase_model.core.species,
    film_phase_model.core.reactions,
    verbose=True,
    check_for_duplicates=True,
)

save_species_dictionary(
    os.path.join(save_directory, "species_dictionary_film.txt"),
    film_phase_model.core.species,
)

print("Getting fragment mapping...")

fragment_based_reaction_mapping = dict()
fragment_based_reaction_label_mapping = dict()
kf_disabled_fragment_based_reactions = []
krev_disabled_fragment_based_reactions = []

label_to_smiles_or_label = dict()
for spc in film_phase_model.core.species:
    label_to_smiles_or_label[spc.label] = spc.smiles + "(L)"
for fragment in fragment_species.values():
    label_to_smiles_or_label[fragment.label] = fragment.label
for oligomer in oligomer_species.values():
    label_to_smiles_or_label[oligomer.label] = oligomer.label
for monomer_label in initial_monomer_labels:
    label_to_smiles_or_label[monomer_label + "(L)"] = monomer_label + "(L)"

# for ind, rxn in enumerate(film_phase_model.core.reactions):
#     label = ""
#     for spc in rxn.reactants:
#         label += f"{spc.label}+"
#     label = label[:-1]
#     label += "<=>"
#     for spc in rxn.products:
#         label += f"{spc.label}+"
#     label = label[:-1]

#     fragment_based_reaction_mapping[label] = dict()

#     if any(spc.label not in fragment_species.keys() and spc.label for spc in rxn.reactants):


for ind, rxn in enumerate(film_phase_model.core.reactions):
    label = ""
    for spc in rxn.reactants:
        label += f"{spc.label}+"
    label = label[:-1]
    label += "<=>"
    for spc in rxn.products:
        label += f"{spc.label}+"
    label = label[:-1]

    fragment_based_reaction_mapping[label] = dict()

    if rxn.family == "R_Addition_MultipleBond":
        krev_disabled_fragment_based_reactions.append(label)
    elif rxn.family == "Diels_alder_addition":
        krev_disabled_fragment_based_reactions.append(label)
    elif rxn.family == "Retroene":
        kf_disabled_fragment_based_reactions.append(label)
    elif rxn.family == "HO2_Elimination_from_PeroxyRadical":
        kf_disabled_fragment_based_reactions.append(label)
    elif rxn.family == "H_Abstraction":
        if any(
            spc.label == "KR" or spc.label == "alkyl_C.(L,oligomer)"
            for spc in rxn.reactants
        ):
            krev_disabled_fragment_based_reactions.append(label)
        elif any(
            spc.label == "KR" or spc.label == "alkyl_C.(L,oligomer)"
            for spc in rxn.products
        ):
            kf_disabled_fragment_based_reactions.append(label)
        # if include_oxygen:
        #     if any(spc.label == "OR" for spc in rxn.reactants):
        #         krev_disabled_fragment_based_reactions.append(label)
        #     elif any(spc.label == "OR" for spc in rxn.products):
        #         kf_disabled_fragment_based_reactions.append(label)
    elif rxn.family == "R_Recombination" and any(
        spc.smiles == "[O][O]" for spc in rxn.reactants
    ):
        krev_disabled_fragment_based_reactions.append(label)

    if rxn.family == "R_Addition_MultipleBond" and any(
        spc.label == "CDB" for spc in rxn.reactants
    ):
        assert len(rxn.products) == 1
        fragment_based_species = rxn.products[0]
        fragment_based_reaction_mapping[label][fragment_based_species.label] = ["KR"]
        fragment_based_reaction_mapping[label][fragment_based_species.label].extend(
            ["CDB"]
            * calculate_subgraph_isomorphisms(
                fragment_based_species.molecule[0], "carbon_double_bond"
            )
        )
        fragment_based_reaction_mapping[label][fragment_based_species.label].extend(
            ["AH"]
            * calculate_subgraph_isomorphisms(
                fragment_based_species.molecule[0], "allylic_CH"
            )
        )
        if include_oxygen:
            fragment_based_reaction_mapping[label][fragment_based_species.label].extend(
                ["HP"]
                * calculate_subgraph_isomorphisms(
                    fragment_based_species.molecule[0], "COOH"
                )
            )
            fragment_based_reaction_mapping[label][fragment_based_species.label].extend(
                ["CP"]
                * calculate_subgraph_isomorphisms(
                    fragment_based_species.molecule[0], "COOC"
                )
            )

        # remove implicitly lost structures
        for spc in rxn.reactants:
            if spc.label in fragment_species_implicit_structures:
                for implicitly_lost_fragment in fragment_species_implicit_structures[
                    spc.label
                ]:
                    if (
                        implicitly_lost_fragment
                        in fragment_based_reaction_mapping[label][
                            fragment_based_species.label
                        ]
                    ):
                        fragment_based_reaction_mapping[label][
                            fragment_based_species.label
                        ].remove(implicitly_lost_fragment)

        fragment_based_reaction_label_mapping[
            label
        ] = f"{label_to_smiles_or_label[rxn.reactants[0].label]} + {label_to_smiles_or_label[rxn.reactants[1].label]} radical addition"

    elif rxn.family == "Diels_alder_addition":
        assert len(rxn.products) == 1
        fragment_based_species = rxn.products[0]
        fragment_based_reaction_mapping[label][fragment_based_species.label] = [
            "CDB"
        ] * calculate_subgraph_isomorphisms(
            fragment_based_species.molecule[0], "carbon_double_bond"
        )
        fragment_based_reaction_mapping[label][fragment_based_species.label].extend(
            ["AH"]
            * calculate_subgraph_isomorphisms(
                fragment_based_species.molecule[0], "allylic_CH"
            )
        )
        if include_oxygen:
            fragment_based_reaction_mapping[label][fragment_based_species.label].extend(
                ["HP"]
                * calculate_subgraph_isomorphisms(
                    fragment_based_species.molecule[0], "COOH"
                )
            )
            fragment_based_reaction_mapping[label][fragment_based_species.label].extend(
                ["CP"]
                * calculate_subgraph_isomorphisms(
                    fragment_based_species.molecule[0], "COOC"
                )
            )

        # remove implicitly lost structures
        for spc in rxn.reactants:
            if spc.label in fragment_species_implicit_structures:
                for implicitly_lost_fragment in fragment_species_implicit_structures[
                    spc.label
                ]:
                    if (
                        implicitly_lost_fragment
                        in fragment_based_reaction_mapping[label][
                            fragment_based_species.label
                        ]
                    ):
                        fragment_based_reaction_mapping[label][
                            fragment_based_species.label
                        ].remove(implicitly_lost_fragment)

        fragment_based_reaction_label_mapping[
            label
        ] = f"{label_to_smiles_or_label[rxn.reactants[0].label]} + {label_to_smiles_or_label[rxn.reactants[1].label]} Diels-Alder addition"

    elif rxn.family == "R_Addition_MultipleBond" and any(
        spc.label in fragment_radical_labels for spc in rxn.reactants
    ):
        assert len(rxn.products) == 1
        if len(rxn.products) == 1:
            fragment_based_species = rxn.products[0]
            if (
                calculate_subgraph_isomorphisms(
                    rxn.products[0].molecule[0], "allylic_C."
                )
                > 0
            ):
                fragment_based_reaction_mapping[label][fragment_based_species.label] = [
                    "AR"
                ]
            else:
                fragment_based_reaction_mapping[label][fragment_based_species.label] = [
                    "KR"
                ]
        elif len(rxn.products) == 2:
            fragment_based_species = rxn.products[1]
            fragment_based_reaction_mapping[label][fragment_based_species.label] = []
        else:
            raise ValueError("Unexpected number of products")

        fragment_based_reaction_mapping[label][fragment_based_species.label].extend(
            ["CD"]
            * calculate_subgraph_isomorphisms(
                fragment_based_species.molecule[0], "conjugated_diene"
            )
        )
        fragment_based_reaction_mapping[label][fragment_based_species.label].extend(
            ["CDB"]
            * (
                calculate_subgraph_isomorphisms(
                    fragment_based_species.molecule[0], "carbon_double_bond"
                )
                - calculate_subgraph_isomorphisms(
                    fragment_based_species.molecule[0], "conjugated_diene"
                )
            )
        )
        fragment_based_reaction_mapping[label][fragment_based_species.label].extend(
            ["AH"]
            * calculate_subgraph_isomorphisms(
                fragment_based_species.molecule[0], "allylic_CH"
            )
        )
        if include_oxygen:
            fragment_based_reaction_mapping[label][fragment_based_species.label].extend(
                ["HP"]
                * calculate_subgraph_isomorphisms(
                    fragment_based_species.molecule[0], "COOH"
                )
            )
            fragment_based_reaction_mapping[label][fragment_based_species.label].extend(
                ["CP"]
                * calculate_subgraph_isomorphisms(
                    fragment_based_species.molecule[0], "COOC"
                )
            )

        # remove implicitly lost structures
        for spc in rxn.reactants:
            if spc.label in fragment_species_implicit_structures:
                for implicitly_lost_fragment in fragment_species_implicit_structures[
                    spc.label
                ]:
                    if (
                        implicitly_lost_fragment
                        in fragment_based_reaction_mapping[label][
                            fragment_based_species.label
                        ]
                    ):
                        fragment_based_reaction_mapping[label][
                            fragment_based_species.label
                        ].remove(implicitly_lost_fragment)

        fragment_based_reaction_label_mapping[
            label
        ] = f"{label_to_smiles_or_label[rxn.reactants[0].label]} + {label_to_smiles_or_label[rxn.reactants[1].label]} radical addition"

    elif rxn.family == "R_Recombination" and any(
        spc.smiles == "[O][O]" for spc in rxn.reactants
    ):
        assert len(rxn.products) == 1
        fragment_based_species = rxn.products[0]
        if fragment_based_species.label != "PR":
            fragment_based_reaction_mapping[label][fragment_based_species.label] = [
                "PR"
            ]
            # fragment_based_reaction_mapping[label][fragment_based_species.label].extend(["CDB"]*calculate_subgraph_isomorphisms(fragment_based_species.molecule[0], "carbon_double_bond"))
            # fragment_based_reaction_mapping[label][fragment_based_species.label].extend(["AH"]*calculate_subgraph_isomorphisms(fragment_based_species.molecule[0], "allylic_CH"))
            # if include_oxygen:
            #     fragment_based_reaction_mapping[label][fragment_based_species.label].extend(["HP"]*calculate_subgraph_isomorphisms(fragment_based_species.molecule[0], "COOH"))
            #     fragment_based_reaction_mapping[label][fragment_based_species.label].extend(["CP"]*calculate_subgraph_isomorphisms(fragment_based_species.molecule[0], "COOC"))

        fragment_based_reaction_label_mapping[
            label
        ] = f"{label_to_smiles_or_label[rxn.reactants[0].label]} + {label_to_smiles_or_label[rxn.reactants[1].label]} <=> PR"

    elif rxn.family == "H_Abstraction":
        fragment_based_reaction_label_mapping[label] = ""
        for spc in rxn.reactants:
            if spc.label in fragment_species or "(L" in spc.label:
                fragment_based_reaction_label_mapping[
                    label
                ] += f"{label_to_smiles_or_label[spc.label]} + "
            else:
                if spc.smiles == "CCCC(C)C":  # alkyl_CH(L,oligomer)
                    inert_label = "inert(L)"
                else:
                    inert_label = "inert(S)"
                    fragment_based_reaction_mapping[label][spc.label] = [inert_label]
                fragment_based_reaction_label_mapping[label] += f"{inert_label} + "
        fragment_based_reaction_label_mapping[
            label
        ] = fragment_based_reaction_label_mapping[label][:-3]
        fragment_based_reaction_label_mapping[label] += " <=> "
        for spc in rxn.products:
            if spc.label in fragment_species or "(L" in spc.label:
                fragment_based_reaction_label_mapping[
                    label
                ] += f"{label_to_smiles_or_label[spc.label]} + "
            else:
                if spc.smiles == "CCCC(C)C":  # alkyl_CH(L,oligomer)
                    inert_label = "inert(L)"
                else:
                    inert_label = "inert(S)"
                    fragment_based_reaction_mapping[label][spc.label] = [inert_label]
                fragment_based_reaction_label_mapping[label] += f"{inert_label} + "
        fragment_based_reaction_label_mapping[
            label
        ] = fragment_based_reaction_label_mapping[label][:-3]
        if (
            "inert(S)" in fragment_based_reaction_label_mapping[label]
            or "inert(L)" in fragment_based_reaction_label_mapping[label]
        ):
            r_label, p_label = fragment_based_reaction_label_mapping[label].split(
                " <=> "
            )
            if "inert(S)" in r_label or "inert(L)" in r_label:
                fragment_based_reaction_label_mapping[label] = f"{p_label} -> {r_label}"
            else:
                fragment_based_reaction_label_mapping[label] = f"{r_label} -> {p_label}"

    elif rxn.family == "Retroene":
        assert len(rxn.reactants) == 1
        fragment_based_species = rxn.reactants[0]
        fragment_based_reaction_mapping[label][fragment_based_species.label] = [
            "CDB"
        ] * calculate_subgraph_isomorphisms(
            fragment_based_species.molecule[0], "carbon_double_bond"
        )
        fragment_based_reaction_mapping[label][fragment_based_species.label].extend(
            ["AH"]
            * calculate_subgraph_isomorphisms(
                fragment_based_species.molecule[0], "allylic_CH"
            )
        )
        if include_oxygen:
            fragment_based_reaction_mapping[label][fragment_based_species.label].extend(
                ["HP"]
                * calculate_subgraph_isomorphisms(
                    fragment_based_species.molecule[0], "COOH"
                )
            )
            fragment_based_reaction_mapping[label][fragment_based_species.label].extend(
                ["CP"]
                * calculate_subgraph_isomorphisms(
                    fragment_based_species.molecule[0], "COOC"
                )
            )

        # remove implicitly lost structures
        for spc in rxn.products:
            if spc.label in fragment_species_implicit_structures:
                for implicitly_lost_fragment in fragment_species_implicit_structures[
                    spc.label
                ]:
                    if (
                        implicitly_lost_fragment
                        in fragment_based_reaction_mapping[label][
                            fragment_based_species.label
                        ]
                    ):
                        fragment_based_reaction_mapping[label][
                            fragment_based_species.label
                        ].remove(implicitly_lost_fragment)

        fragment_based_reaction_label_mapping[
            label
        ] = f"{label_to_smiles_or_label[rxn.products[0].label]} + {label_to_smiles_or_label[rxn.products[1].label]} ene reaction"

    elif rxn.family == "HO2_Elimination_from_PeroxyRadical":
        assert len(rxn.reactants) == 1
        fragment_based_species = rxn.reactants[0]
        fragment_based_reaction_mapping[label][fragment_based_species.label] = ["PR"]

        fragment_based_reaction_label_mapping[label] = "HO2 + CDB concerted addition"
    else:
        fragment_based_reaction_label_mapping[label] = label

    print(rxn)
    print(fragment_based_reaction_mapping[label])
    print(f"kf disabled: {label in kf_disabled_fragment_based_reactions}")
    print(f"krev disabled: {label in krev_disabled_fragment_based_reactions}")

path = os.path.join(save_directory, "chem_film_phase.rms")
with open(path, "r") as f:
    film_mech_dict = yaml.load(f, Loader=yaml.SafeLoader)

fragment_intermediate = []
for rxn_dict in film_mech_dict["Reactions"]:
    rxn_string = ""
    for spc in rxn_dict["reactants"]:
        rxn_string += f"{spc}+"
    rxn_string = rxn_string[:-1]
    rxn_string += "<=>"
    for spc in rxn_dict["products"]:
        rxn_string += f"{spc}+"
    rxn_string = rxn_string[:-1]
    if rxn_string in kf_disabled_fragment_based_reactions:
        rxn_dict["forwardable"] = False
    if rxn_string in krev_disabled_fragment_based_reactions:
        rxn_dict["reversible"] = False
    if rxn_string in fragment_based_reaction_mapping:
        mapping = fragment_based_reaction_mapping[rxn_string]
        rxn_dict["fragmentbasedreactants"] = []
        rxn_dict["fragmentbasedproducts"] = []
        for spc in rxn_dict["reactants"]:
            if spc in mapping:
                fragment_intermediate.append(spc)
                rxn_dict["fragmentbasedreactants"].extend(mapping[spc])
            else:
                rxn_dict["fragmentbasedreactants"].append(spc)
        for spc in rxn_dict["products"]:
            if spc in mapping:
                fragment_intermediate.append(spc)
                rxn_dict["fragmentbasedproducts"].extend(mapping[spc])
            else:
                rxn_dict["fragmentbasedproducts"].append(spc)
    if rxn_string in fragment_based_reaction_label_mapping:
        rxn_dict["comment"] = fragment_based_reaction_label_mapping[rxn_string]

fragment_intermediate = list(set(fragment_intermediate))
for spc_dict in film_mech_dict["Phases"][0]["Species"]:
    if spc_dict["name"] in fragment_species:
        spc_dict["isfragment"] = True
    if spc_dict["name"] in fragment_intermediate:
        spc_dict["isfragmentintermediate"] = True

film_spc_dicts = [
    spc_dict
    for spc_dict in film_mech_dict["Phases"][0]["Species"]
    if spc_dict.get("isfragmentintermediate", False)
    or spc_dict.get("isfragment", False)
]
liq_spc_dicts = [
    spc_dict
    for spc_dict in film_mech_dict["Phases"][0]["Species"]
    if not spc_dict.get("isfragmentintermediate", False)
    and not spc_dict.get("isfragment", False)
]

film_mech_dict["Phases"][0]["Species"] = film_spc_dicts
film_mech_dict["Phases"][0]["name"] = "film"
film_mech_dict["Phases"].append({})
film_mech_dict["Phases"][1]["Species"] = liq_spc_dicts
film_mech_dict["Phases"][1]["name"] = "liquid"

path = os.path.join(save_directory, "chem_film_phase.rms")
with open(path, "w+") as f:
    yaml.dump(film_mech_dict, stream=f)


print("Getting species labels relevant to ASF distirbution and initial conditions...")

monomer_labels = initial_monomer_labels[:]
monomer_allylic_C_labels = []
monomer_alkyl_C_labels = []

if include_oxygen:
    monomer_COO_labels = []

for rxn in liqrxns:
    if rxn.family == "Disproportionation":
        if all(spc.label in initial_monomer_labels for spc in rxn.products):
            for spc in rxn.reactants:
                if spc.multiplicity == 2:
                    if (
                        calculate_subgraph_isomorphisms(spc.molecule[0], "allylic_C.")
                        > 0
                    ):
                        monomer_allylic_C_labels.append(spc.label)
                    elif (
                        calculate_subgraph_isomorphisms(spc.molecule[0], "alkyl_C.") > 0
                    ):
                        monomer_alkyl_C_labels.append(spc.label)
                    else:
                        continue
                    mol = spc.molecule[0].copy(deep=True)
                    mol.saturate_radicals()
                    if mol.smiles in liq_smiles_to_spc_dict:
                        label = liq_smiles_to_spc_dict[mol.smiles].label
                        monomer_labels.append(label)
monomer_labels = list(set(monomer_labels))
monomer_allylic_C_labels = list(set(monomer_allylic_C_labels))
monomer_alkyl_C_labels = list(set(monomer_alkyl_C_labels))

if include_oxygen:
    monomer_carbon_radical_labels = monomer_allylic_C_labels + monomer_alkyl_C_labels
    for rxn in liqrxns:
        if rxn.family == "R_Recombination":
            if any(spc.smiles == "[O][O]" for spc in rxn.reactants) and any(
                spc.label in monomer_carbon_radical_labels for spc in rxn.reactants
            ):
                monomer_COO_labels.extend([spc.label for spc in rxn.products])

    monomer_COO_labels = list(set(monomer_COO_labels))

print("Monomer smiles...")
for label in monomer_labels:
    print(label)
    spc = liq_label_to_spc_dict[label]
    print(spc.smiles)

print("Monomer allylic radical smiles...")
for label in monomer_allylic_C_labels:
    print(label)
    spc = liq_label_to_spc_dict[label]
    print(spc.smiles)

print("Monomer alkyl radical smiles...")
for label in monomer_alkyl_C_labels:
    print(label)
    spc = liq_label_to_spc_dict[label]
    print(spc.smiles)

if include_oxygen:
    print("Monomer peroxyl radical smiles...")
    for label in monomer_COO_labels:
        print(label)
        spc = liq_label_to_spc_dict[label]
        print(spc.smiles)

dimer_labels = []
dimer_allylic_C_labels = []
dimer_alkyl_C_labels = []
if include_oxygen:
    dimer_COO_labels = []
    dimer_COOH_labels = []

for rxn in liqrxns:
    if rxn.family == "Diels_alder_addition":
        if all([spc.label in monomer_labels for spc in rxn.reactants]):
            dimer_labels.append(rxn.products[0].label)
    elif rxn.family == "Retroene":
        if all([spc.label in monomer_labels for spc in rxn.products]):
            dimer_labels.append(rxn.reactants[0].label)
    elif rxn.family == "R_Recombination":
        if all(
            spc.label in monomer_allylic_C_labels or spc.label in monomer_alkyl_C_labels
            for spc in rxn.reactants
        ):
            dimer_labels.append(rxn.products[0].label)
    elif rxn.family == "R_Addition_MultipleBond":
        if len(rxn.reactants) == 2:
            two_spc_side = rxn.reactants
            one_spc_side = rxn.products
        elif len(rxn.products) == 2:
            two_spc_side = rxn.products
            one_spc_side = rxn.reactants
        else:
            continue
        if any(
            [
                spc.label in monomer_allylic_C_labels
                or spc.label in monomer_alkyl_C_labels
                for spc in two_spc_side
            ]
        ) and any([spc.label in monomer_labels for spc in two_spc_side]):
            if one_spc_side[0].multiplicity == 2:
                if (
                    calculate_subgraph_isomorphisms(
                        one_spc_side[0].molecule[0], "allylic_C."
                    )
                    > 0
                ):
                    dimer_allylic_C_labels.append(one_spc_side[0].label)
                elif (
                    calculate_subgraph_isomorphisms(
                        one_spc_side[0].molecule[0], "alkyl_C."
                    )
                    > 0
                ):
                    dimer_alkyl_C_labels.append(one_spc_side[0].label)
                else:
                    continue
                mol = one_spc_side[0].molecule[0].copy(deep=True)
                mol.saturate_radicals()
                if mol.smiles in liq_smiles_to_spc_dict:
                    dimer_labels.append(liq_smiles_to_spc_dict[mol.smiles].label)

dimer_labels = list(set(dimer_labels))
dimer_allylic_C_labels = list(set(dimer_allylic_C_labels))
dimer_alkyl_C_labels = list(set(dimer_alkyl_C_labels))

if include_oxygen:
    dimer_carbon_radical_labels = dimer_allylic_C_labels + dimer_alkyl_C_labels
    for rxn in liqrxns:
        if rxn.family == "R_Recombination":
            if any(spc.smiles == "[O][O]" for spc in rxn.reactants) and any(
                spc.label in dimer_carbon_radical_labels for spc in rxn.reactants
            ):
                dimer_COO_labels.extend([spc.label for spc in rxn.products])
                for spc in rxn.products:
                    mol = spc.molecule[0].copy(deep=True)
                    mol.saturate_radicals()
                    if mol.smiles in liq_smiles_to_spc_dict:
                        label = liq_smiles_to_spc_dict[mol.smiles].label
                    dimer_COOH_labels.append(label)

    dimer_COO_labels = list(set(dimer_COO_labels))
    dimer_COOH_labels = list(set(dimer_COOH_labels))

print("Dimer smiles...")
for label in dimer_labels:
    spc = liq_label_to_spc_dict[label]
    print(spc.smiles)

print("Dimer allylic radical smiles...")
for label in dimer_allylic_C_labels:
    spc = liq_label_to_spc_dict[label]
    print(spc.smiles)

print("Dimer alkyl radical smiles...")
for label in dimer_alkyl_C_labels:
    spc = liq_label_to_spc_dict[label]
    print(spc.smiles)

if include_oxygen:
    print("Dimer COO smiles...")
    for label in dimer_COO_labels:
        spc = liq_label_to_spc_dict[label]
        print(spc.smiles)

    print("Dimer COOH smiles...")
    for label in dimer_COOH_labels:
        spc = liq_label_to_spc_dict[label]
        print(spc.smiles)

liquid_species_mapping = dict()
liquid_species_mapping["dimer(L)"] = dimer_labels
liquid_species_mapping["dimer_allylic_C.(L)"] = dimer_allylic_C_labels
liquid_species_mapping["dimer_alkyl_C.(L)"] = dimer_alkyl_C_labels
if include_oxygen:
    liquid_species_mapping["dimer_COO.(L)"] = dimer_COO_labels
    liquid_species_mapping["dimer_COOH(L)"] = dimer_COOH_labels
liquid_species_mapping["monomer(L)"] = monomer_labels

liquid_species_mapping["allylic_CH(L)"] = list()
liquid_species_mapping["allylic_C.(L)"] = list()
liquid_species_mapping["alkyl_C.(L)"] = list()
liquid_species_mapping["C=C(L)"] = list()
liquid_species_mapping["conjugated_diene(L)"] = list()
if include_oxygen:
    liquid_species_mapping["COO.(L)"] = list()
    liquid_species_mapping["COOC(L)"] = list()
    liquid_species_mapping["COOH(L)"] = list()
    liquid_species_mapping["CO.(L)"] = list()

for spc in liqspcs:
    liquid_species_mapping["allylic_CH(L)"].extend(
        [spc.label] * calculate_subgraph_isomorphisms(spc.molecule[0], "allylic_CH")
    )
    liquid_species_mapping["C=C(L)"].extend(
        [spc.label]
        * calculate_subgraph_isomorphisms(spc.molecule[0], "carbon_double_bond")
    )

    if spc.multiplicity == 2:
        liquid_species_mapping["allylic_C.(L)"].extend(
            [spc.label] * calculate_subgraph_isomorphisms(spc.molecule[0], "allylic_C.")
        )
        liquid_species_mapping["alkyl_C.(L)"].extend(
            [spc.label] * calculate_subgraph_isomorphisms(spc.molecule[0], "alkyl_C.")
        )
        # film conjugated dienes come from radicals
        liquid_species_mapping["conjugated_diene(L)"].extend(
            [spc.label]
            * calculate_subgraph_isomorphisms(spc.molecule[0], "conjugated_diene")
        )

    if include_oxygen:
        if spc.multiplicity == 2:
            liquid_species_mapping["COO.(L)"].extend(
                [spc.label] * calculate_subgraph_isomorphisms(spc.molecule[0], "COO.")
            )
            liquid_species_mapping["CO.(L)"].extend(
                [spc.label] * calculate_subgraph_isomorphisms(spc.molecule[0], "CO.")
            )
        liquid_species_mapping["COOC(L)"].extend(
            [spc.label] * calculate_subgraph_isomorphisms(spc.molecule[0], "COOC")
        )
        liquid_species_mapping["COOH(L)"].extend(
            [spc.label] * calculate_subgraph_isomorphisms(spc.molecule[0], "COOH")
        )

path = os.path.join(save_directory, f"liquid_species_mapping.yml")
with open(path, "w+") as f:
    yaml.dump(liquid_species_mapping, stream=f)

end = time.time()
print(f"Time taken: {end-start} seconds")
print("Done!")
