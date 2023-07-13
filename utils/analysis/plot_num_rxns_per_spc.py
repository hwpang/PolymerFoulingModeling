import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed

from rmgpy.chemkin import load_chemkin_file
from rmgpy import settings
from rmgpy.rmg.main import RMG

# Load liquid-phase model
debutanizer_liq_chemkin_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/liquid_mechanism/chem_annotated.inp"
debutanizer_liq_spc_dict_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/liquid_mechanism/species_dictionary.txt"
debutanizer_film_chemkin_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/film_mechanism/chem_annotated_film.inp"
debutanizer_film_spc_dict_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/film_mechanism/species_dictionary_film.txt"

QCMD_liq_chemkin_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/QCMD_cell_model/liquid_mechanism/chem_annotated.inp"
QCMD_liq_spc_dict_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/QCMD_cell_model/liquid_mechanism/species_dictionary.txt"
QCMD_film_chemkin_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/QCMD_cell_model/film_mechanism/chem_annotated_film.inp"
QCMD_film_spc_dict_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/QCMD_cell_model/film_mechanism/species_dictionary_film.txt"

def plot_num_rxns_per_spc(spcs, rxns, ax=None):
    """
    Plot the number of reactions per species
    """
    if ax is None:
        fig, ax = plt.subplots()
    num_rxns_per_spc = np.zeros(len(spcs))
    for i, spc in enumerate(spcs):
        for rxn in rxns:
            if spc in rxn.reactants or spc in rxn.products:
                num_rxns_per_spc[i] += 1
    ax.hist(num_rxns_per_spc,
            bins=20,
            edgecolor="black",)

    inds = np.argsort(num_rxns_per_spc)[::-1]
    for ind in inds[:5]:
        print("Species with the most reactions: {}".format(spcs[ind].label))

def plot_num_rxns_per_family(rxns, ax=None):
    """
    Plot the number of reactions per family
    """
    if ax is None:
        fig, ax = plt.subplots()
    families = []
    for rxn in rxns:
        if rxn.family not in families:
            families.append(rxn.family)
    num_rxns_per_family = np.zeros(len(families))
    for i, family in enumerate(families):
        for rxn in rxns:
            if rxn.family == family:
                num_rxns_per_family[i] += 1
    # only plot the top N families and lump the rest into "other"
    N = 5
    inds = np.argsort(num_rxns_per_family)[::-1][:N]
    families_sorted = [families[ind] for ind in inds] + ["other"]
    xs = np.arange(len(families_sorted))
    ys = list(num_rxns_per_family[inds]) + [np.sum(num_rxns_per_family) - np.sum(num_rxns_per_family[inds])]
    print(xs)
    print(ys)
    ax.bar(xs, ys)
    ax.set_xticks(xs)
    ax.set_xticklabels(families_sorted, rotation=45, ha="right")

def plot_num_rxns_per_spc_by_spc_size(spcs, rxns, ax=None):
    """
    Plot the number of reactions per species, grouped by species size
    """
    if ax is None:
        fig, ax = plt.subplots()
    num_rxns_per_spc = np.zeros(len(spcs))
    for i, spc in enumerate(spcs):
        for rxn in rxns:
            if spc in rxn.reactants or spc in rxn.products:
                num_rxns_per_spc[i] += 1
    num_rxns_per_spc_grouped_by_spc_size = {}
    for i, spc in enumerate(spcs):
        num_heavy_atoms = spc.molecule[0].get_num_atoms() - spc.molecule[0].get_num_atoms('H')
        if num_heavy_atoms not in num_rxns_per_spc_grouped_by_spc_size:
            num_rxns_per_spc_grouped_by_spc_size[num_heavy_atoms] = []
        num_rxns_per_spc_grouped_by_spc_size[num_heavy_atoms].append(num_rxns_per_spc[i])
    # plot violin plot
    ax.violinplot(num_rxns_per_spc_grouped_by_spc_size.values(),
                  num_rxns_per_spc_grouped_by_spc_size.keys())


fragments = set(["CDB", "AR", "KR", "AH", "PR", "CP", "OR", "OH"])

debutanizer_liqspcs, debutanizer_liqrxns = load_chemkin_file(debutanizer_liq_chemkin_path, debutanizer_liq_spc_dict_path)
print(f"Debutanizer liquid mechanism: {len(debutanizer_liqspcs)} species, {len(debutanizer_liqrxns)} reactions")
debutanizer_filmspcs, debutanizer_filmrxns = load_chemkin_file(debutanizer_film_chemkin_path, debutanizer_film_spc_dict_path)
debutanizer_filmspcs = [spc for spc in debutanizer_filmspcs if spc.label in fragments or "(L)" in spc.label]
print(f"Debutanizer film mechanism: {len(debutanizer_filmspcs)} species, {len(debutanizer_filmrxns)} reactions")

QCMD_liqspcs, QCMD_liqrxns = load_chemkin_file(QCMD_liq_chemkin_path, QCMD_liq_spc_dict_path)
print(f"QCMD liquid mechanism: {len(QCMD_liqspcs)} species, {len(QCMD_liqrxns)} reactions")
QCMD_filmspcs, QCMD_filmrxns = load_chemkin_file(QCMD_film_chemkin_path, QCMD_film_spc_dict_path)
QCMD_filmspcs = [spc for spc in QCMD_filmspcs if spc.label in fragments or "(L)" in spc.label]
print(f"QCMD film mechanism: {len(QCMD_filmspcs)} species, {len(QCMD_filmrxns)} reactions")

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

for rxn in debutanizer_liqrxns+debutanizer_filmrxns+QCMD_liqrxns+QCMD_filmrxns:
    rxn.family = get_reaction_family(rxn)

fig, axs = plt.subplots(2, 2, figsize=(8, 8))

ax = axs[0, 0]
plot_num_rxns_per_family(debutanizer_liqrxns, ax=ax)
ax.set_title("Debutanizer liquid submodel")
ax.set_ylabel("Count")

ax = axs[0, 1]
plot_num_rxns_per_family(debutanizer_filmrxns, ax=ax)
ax.set_title("Debutanizer film submodel")

ax = axs[1, 0]
plot_num_rxns_per_family(QCMD_liqrxns, ax=ax)
ax.set_title("QCMD liquid submodel")
ax.set_ylabel("Count")
ax.set_xlabel("Reaction family")

ax = axs[1, 1]
plot_num_rxns_per_family(QCMD_filmrxns, ax=ax)
ax.set_title("QCMD film submodel")
ax.set_xlabel("Reaction family")

fig.tight_layout()
fig.savefig("Figures/num_rxns_per_family.pdf", bbox_inches="tight")

fig, axs = plt.subplots(2, 2, figsize=(8, 6))

ax = axs[0, 0]
plot_num_rxns_per_spc(debutanizer_liqspcs, debutanizer_liqrxns, ax=ax)
ax.set_title("Debutanizer liquid submodel")
ax.set_ylabel("Number of species")

ax = axs[0, 1]
plot_num_rxns_per_spc(debutanizer_filmspcs, debutanizer_filmrxns, ax=ax)
ax.set_title("Debutanizer film submodel")

ax = axs[1, 0]
plot_num_rxns_per_spc(QCMD_liqspcs, QCMD_liqrxns, ax=ax)
ax.set_title("QCMD liquid submodel")
ax.set_ylabel("Number of species")
ax.set_xlabel("Number of involved reactions")

ax = axs[1, 1]
plot_num_rxns_per_spc(QCMD_filmspcs, QCMD_filmrxns, ax=ax)
ax.set_title("QCMD film submodel")
ax.set_xlabel("Number of involved reactions")

fig.tight_layout()
fig.savefig("Figures/num_rxns_per_spc.pdf", bbox_inches="tight")
