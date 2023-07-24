import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed

from rmgpy.chemkin import load_chemkin_file
from rmgpy import settings
from rmgpy.rmg.main import RMG
from rmgpy.data.kinetics.library import LibraryReaction

# Load liquid-phase model
debutanizer_liq_chemkin_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/liquid_mechanism/chem_annotated.inp"
debutanizer_liq_spc_dict_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/liquid_mechanism/species_dictionary.txt"
debutanizer_film_chemkin_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/film_mechanism/chem_annotated_film.inp"
debutanizer_film_spc_dict_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/film_mechanism/species_dictionary_film.txt"

QCMD_liq_chemkin_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/QCMD_cell_model/liquid_mechanism/chem_annotated.inp"
QCMD_liq_spc_dict_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/QCMD_cell_model/liquid_mechanism/species_dictionary.txt"
QCMD_film_chemkin_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/QCMD_cell_model/film_mechanism/chem_annotated_film.inp"
QCMD_film_spc_dict_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/QCMD_cell_model/film_mechanism/species_dictionary_film.txt"

def categorize_spc_by_thermo_source(spcs, ax=None):
    spc_by_thermo_source = {}
    spc_by_thermo_source["Library"] = []
    spc_by_thermo_source["QM"] = []
    spc_by_thermo_source["GAV"] = []

    QM_libraries = set(['hwpang_fouling',
        'thermo_combined_new',
        'Conjugated_diene',
        'ExptModelTraceO2',
        'multi_trays_combined',
        'multi_trays_v3',
        'multi_trays',
        'thermo_combined',])

    for spc in spcs:
        if any(QM_library in spc.thermo.comment for QM_library in QM_libraries):
            spc_by_thermo_source["QM"].append(spc)
        elif "Thermo library" in spc.thermo.comment:
            spc_by_thermo_source["Library"].append(spc)
        elif "Thermo group additivity estimation" in spc.thermo.comment:
            spc_by_thermo_source["GAV"].append(spc)
        else:
            print("Unknown thermo source for species {0}".format(spc.label))

    print("Species thermo source:")
    print(f"Total: {len(spcs)}")
    for source, spc_list in spc_by_thermo_source.items():
        count = len(spc_list)
        print(f"{source}: {count} ({count/len(spcs)*100:.0f}%)")

    GAV_spc_sizes = [spc.molecule[0].get_num_atoms() - spc.molecule[0].get_num_atoms('H') for spc in spc_by_thermo_source["GAV"]]
    QM_spc_sizes = [spc.molecule[0].get_num_atoms() - spc.molecule[0].get_num_atoms('H') for spc in spc_by_thermo_source["QM"]]
    Library_spc_sizes = [spc.molecule[0].get_num_atoms() - spc.molecule[0].get_num_atoms('H') for spc in spc_by_thermo_source["Library"]]

    ax.hist(
        [GAV_spc_sizes, QM_spc_sizes, Library_spc_sizes],
        bins=np.array(range(25))-0.5,
        label=[f"GAV (n={len(GAV_spc_sizes)})", f"QM (n={len(QM_spc_sizes)})", f"Library (n={len(Library_spc_sizes)})"],
        stacked=True,
        edgecolor='black',
    )
    ax.legend(bbox_to_anchor=(0., -0.1, 1., -0.1), loc=1)

def categorize_rxn_by_kinetic_source(rxns):
    rxn_by_kinetic_source = {}
    rxn_by_kinetic_source["Library"] = []
    rxn_by_kinetic_source["QM"] = []
    rxn_by_kinetic_source["Rate rules"] = []

    QM_libraries = set(['hwpang_fouling','Conjugated_diene'])

    for rxn in rxns:
        if isinstance(rxn, LibraryReaction):
            if "Fitted to " in rxn.kinetics.comment:
                rxn_by_kinetic_source["QM"].append(rxn)
            else:
                rxn_by_kinetic_source["Library"].append(rxn)
        else:
            if "From training reaction" in rxn.kinetics.comment or "Matched reaction" in rxn.kinetics.comment:
                rxn_by_kinetic_source["Library"].append(rxn)
            elif "Estimated" in rxn.kinetics.comment:
                rxn_by_kinetic_source["Rate rules"].append(rxn)
            else:
                rxn_by_kinetic_source["Rate rules"].append(rxn) 

        # if any(QM_library in rxn.kinetics.comment for QM_library in QM_libraries) or "data points" in rxn.kinetics.comment:
        #     rxn_by_kinetic_source["QM"].append(rxn)
        # elif "From training reaction" in rxn.kinetics.comment or "Matched reaction" in rxn.kinetics.comment:
        #     rxn_by_kinetic_source["Library"].append(rxn)
        # elif "Estimated" in rxn.kinetics.comment:
        #     rxn_by_kinetic_source["Rate rules"].append(rxn)
        # else:
        #     print(f"Unknown kinetic source for reaction {rxn.index} {rxn}, rxn.kinetics.comment: {rxn.kinetics.comment}, rxn.family: {rxn.family}, type(rxn): {type(rxn)}")
        #     rxn_by_kinetic_source["Rate rules"].append(rxn)

    
    print("Reaction kinetic source:")
    print(f"Total: {len(rxns)}")
    for source, rxn_list in rxn_by_kinetic_source.items():
        count = len(rxn_list)
        print(f"{source}: {count} ({count/len(rxns)*100:.0f}%)")

    return rxn_by_kinetic_source

def categorize_rxn_by_kinetic_source_and_size(rxns, ax=None):

    rxn_by_kinetic_source = categorize_rxn_by_kinetic_source(rxns)
    
    GAV_rxn_sizes = [sum(spc.molecule[0].get_num_atoms() - spc.molecule[0].get_num_atoms('H') for spc in rxn.reactants) for rxn in rxn_by_kinetic_source["Rate rules"]]
    QM_rxn_sizes = [sum(spc.molecule[0].get_num_atoms() - spc.molecule[0].get_num_atoms('H') for spc in rxn.reactants) for rxn in rxn_by_kinetic_source["QM"]]
    Library_rxn_sizes = [sum(spc.molecule[0].get_num_atoms() - spc.molecule[0].get_num_atoms('H') for spc in rxn.reactants) for rxn in rxn_by_kinetic_source["Library"]]

    ax.hist(
        [GAV_rxn_sizes, QM_rxn_sizes, Library_rxn_sizes],
        bins=np.array(range(35))-0.5,
        label=[f"Rate rules (n={len(GAV_rxn_sizes)})", f"QM (n={len(QM_rxn_sizes)})", f"Library/Training (n={len(Library_rxn_sizes)})"],
        stacked=True,
        edgecolor='black',
    )
    ax.legend(bbox_to_anchor=(0., -0.1, 1., -0.1), loc=1)

def categprize_rxn_by_family(rxns):
    families = set(rxn.family for rxn in rxns)
    rxn_by_family = {family: [] for family in families}
    for rxn in rxns:
        rxn_by_family[rxn.family].append(rxn)
    return rxn_by_family

def categorize_rxn_by_kinetic_source_and_family(rxns, ax=None):
    rxn_by_family = categprize_rxn_by_family(rxns)
    families = list(rxn_by_family.keys())
    num_rxns_by_family = [len(rxn_by_family[family]) for family in families]

    N = 5
    inds = np.argsort(num_rxns_by_family)[::-1][:N]
    families_sorted = [families[ind] for ind in inds] + ["other"]

    other_rxns = []
    for family in list(rxn_by_family.keys()):
        if family not in families_sorted:
            other_rxns += rxn_by_family[family]
            rxn_by_family.pop(family)
    rxn_by_family["other"] = other_rxns

    rxn_by_kinetic_source_and_family = {}
    for family, rxn_list in rxn_by_family.items():
        rxn_by_kinetic_source_and_family[family] = categorize_rxn_by_kinetic_source(rxn_list)
    
    xs = np.arange(len(families_sorted))
    bottom = np.zeros(len(families_sorted))
    ys = np.array([len(rxn_by_kinetic_source_and_family[family]["Rate rules"]) for family in families_sorted])
    ax.bar(xs, ys, bottom=bottom, label="Rate rules", edgecolor='black')
    bottom += ys

    ys = np.array([len(rxn_by_kinetic_source_and_family[family]["QM"]) for family in families_sorted])
    ax.bar(xs, ys, bottom=bottom, label="QM", edgecolor='black')
    bottom += ys

    ys = np.array([len(rxn_by_kinetic_source_and_family[family]["Library"]) for family in families_sorted])
    ax.bar(xs, ys, bottom=bottom, label="Library/Training", edgecolor='black')

    ax.set_xticks(xs)
    ax.set_xticklabels(families_sorted, rotation=45, ha="right")

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

fragments = set(["CDB", "AR", "KR", "AH", "PR", "CP", "OR", "OH"])

fig, axs = plt.subplots(3, 2, figsize=(8, 11))

liqspcs, liqrxns = load_chemkin_file(debutanizer_liq_chemkin_path, debutanizer_liq_spc_dict_path)
filmspcs, filmrxns = load_chemkin_file(debutanizer_film_chemkin_path, debutanizer_film_spc_dict_path)

allspcs = liqspcs + [spc for spc in filmspcs if spc.label in fragments]
allrxns = liqrxns + filmrxns

ax = axs[0, 0]
categorize_spc_by_thermo_source(allspcs, ax=ax)
ax.set_title("(a) Debutanizer model", loc="left")
ax.set_ylabel("Number of species")
ax.set_xlabel("Number of heavy atoms")

ax = axs[1, 0]
categorize_rxn_by_kinetic_source_and_size(allrxns, ax=ax)
ax.set_title("(b) Debutanizer model", loc="left")
ax.set_ylabel("Number of reactions")
ax.set_xlabel("Number of heavy atoms")

for rxn in allrxns:
    rxn.family = get_reaction_family(rxn)

ax = axs[2, 0]
categorize_rxn_by_kinetic_source_and_family(allrxns, ax=ax)
ax.set_title("(c) Debutanizer model", loc="left")
ax.set_ylabel("Number of reactions")
ax.set_xlabel("Reaction family")

liqspcs, liqrxns = load_chemkin_file(QCMD_liq_chemkin_path, QCMD_liq_spc_dict_path)
filmspcs, filmrxns = load_chemkin_file(QCMD_film_chemkin_path, QCMD_film_spc_dict_path)

allspcs = liqspcs + [spc for spc in filmspcs if spc.label in fragments]
allrxns = liqrxns + filmrxns

ax = axs[0, 1]
categorize_spc_by_thermo_source(allspcs, ax=ax)
ax.set_title("(c) QCMD cell model", loc="left")
ax.set_xlabel("Number of heavy atoms")

ax = axs[1, 1]
categorize_rxn_by_kinetic_source_and_size(allrxns, ax=ax)
ax.set_title("(d) QCMD cell model", loc="left")
ax.set_xlabel("Number of heavy atoms")

for rxn in allrxns:
    rxn.family = get_reaction_family(rxn)

ax = axs[2, 1]
categorize_rxn_by_kinetic_source_and_family(allrxns, ax=ax)
ax.set_title("(e) QCMD cell model", loc="left")
ax.set_xlabel("Reaction family")

fig.tight_layout()
fig.savefig("Figures/thermo_and_kinetic_source.pdf", bbox_inches='tight')


