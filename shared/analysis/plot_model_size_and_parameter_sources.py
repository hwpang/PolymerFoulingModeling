import argparse
import numpy as np
import matplotlib.pyplot as plt

from rmgpy.chemkin import load_chemkin_file
from rmgpy import settings
from rmgpy.rmg.main import RMG
from rmgpy.data.kinetics import LibraryReaction

from utils import get_reaction_family


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model_name",
        type=str,
        required=True,
        help="The name of the model.",
    )

    args = parser.parse_args()
    model_name = args.model_name

    return (model_name,)


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
    ax.hist(
        num_rxns_per_spc,
        bins=20,
        edgecolor="black",
    )

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
    ys = list(num_rxns_per_family[inds]) + [
        np.sum(num_rxns_per_family) - np.sum(num_rxns_per_family[inds])
    ]
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
        num_heavy_atoms = spc.molecule[0].get_num_atoms() - spc.molecule[
            0
        ].get_num_atoms("H")
        if num_heavy_atoms not in num_rxns_per_spc_grouped_by_spc_size:
            num_rxns_per_spc_grouped_by_spc_size[num_heavy_atoms] = []
        num_rxns_per_spc_grouped_by_spc_size[num_heavy_atoms].append(
            num_rxns_per_spc[i]
        )
    # plot violin plot
    ax.violinplot(
        num_rxns_per_spc_grouped_by_spc_size.values(),
        num_rxns_per_spc_grouped_by_spc_size.keys(),
    )


def categorize_spc_by_thermo_source(spcs):
    spc_by_thermo_source = {}
    spc_by_thermo_source["Library"] = []
    spc_by_thermo_source["QM"] = []
    spc_by_thermo_source["GAV"] = []

    QM_libraries = set(
        [
            "hwpang_fouling",
            "thermo_combined_new",
            "Conjugated_diene",
            "ExptModelTraceO2",
            "multi_trays_combined",
            "multi_trays_v3",
            "multi_trays",
            "thermo_combined",
        ]
    )

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

    return spc_by_thermo_source

def get_num_heavy_atoms(spc, use_carbon=False, use_oxygen=False):
    assert not (use_carbon and use_oxygen)

    if use_carbon:
        return spc.molecule[0].get_num_atoms("C")
    elif use_oxygen:
        return spc.molecule[0].get_num_atoms("O")
    else:
        return spc.molecule[0].get_num_atoms() - spc.molecule[0].get_num_atoms("H")

def plot_spc_thermo_source_by_spc_size(spc_by_thermo_source, ax=None, use_carbon=False, use_oxygen=False):

    GAV_spc_sizes = [
        get_num_heavy_atoms(spc, use_carbon=use_carbon, use_oxygen=use_oxygen)
        for spc in spc_by_thermo_source["GAV"]
    ]
    QM_spc_sizes = [
        get_num_heavy_atoms(spc, use_carbon=use_carbon, use_oxygen=use_oxygen)
        for spc in spc_by_thermo_source["QM"]
    ]
    Library_spc_sizes = [
        get_num_heavy_atoms(spc, use_carbon=use_carbon, use_oxygen=use_oxygen)
        for spc in spc_by_thermo_source["Library"]
    ]

    ax.hist(
        [GAV_spc_sizes, QM_spc_sizes, Library_spc_sizes],
        bins=np.array(range(25)) - 0.5,
        label=[
            f"GAV (n={len(GAV_spc_sizes)})",
            f"QM (n={len(QM_spc_sizes)})",
            f"Library (n={len(Library_spc_sizes)})",
        ],
        stacked=True,
        edgecolor="black",
    )


def categorize_rxn_by_kinetic_source(rxns):
    rxn_by_kinetic_source = {}
    rxn_by_kinetic_source["Library"] = []
    rxn_by_kinetic_source["QM"] = []
    rxn_by_kinetic_source["Rate rules"] = []

    QM_libraries = set(["hwpang_fouling", "Conjugated_diene"])

    for rxn in rxns:
        if isinstance(rxn, LibraryReaction):
            if "Fitted to " in rxn.kinetics.comment:
                rxn_by_kinetic_source["QM"].append(rxn)
            else:
                rxn_by_kinetic_source["Library"].append(rxn)
        else:
            if (
                "From training reaction" in rxn.kinetics.comment
                or "Matched reaction" in rxn.kinetics.comment
            ):
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

    GAV_rxn_sizes = [
        sum(
            get_num_heavy_atoms(spc)
            for spc in rxn.reactants
        )
        for rxn in rxn_by_kinetic_source["Rate rules"]
    ]
    QM_rxn_sizes = [
        sum(
            get_num_heavy_atoms(spc)
            for spc in rxn.reactants
        )
        for rxn in rxn_by_kinetic_source["QM"]
    ]
    Library_rxn_sizes = [
        sum(
            get_num_heavy_atoms(spc)
            for spc in rxn.reactants
        )
        for rxn in rxn_by_kinetic_source["Library"]
    ]

    ax.hist(
        [GAV_rxn_sizes, QM_rxn_sizes, Library_rxn_sizes],
        bins=np.array(range(35)) - 0.5,
        label=[
            f"Rate rules (n={len(GAV_rxn_sizes)})",
            f"QM (n={len(QM_rxn_sizes)})",
            f"Library/Training (n={len(Library_rxn_sizes)})",
        ],
        stacked=True,
        edgecolor="black",
    )
    ax.legend(bbox_to_anchor=(0.0, -0.1, 1.0, -0.1), loc=1)


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
        rxn_by_kinetic_source_and_family[family] = categorize_rxn_by_kinetic_source(
            rxn_list
        )

    xs = np.arange(len(families_sorted))
    bottom = np.zeros(len(families_sorted))
    ys = np.array(
        [
            len(rxn_by_kinetic_source_and_family[family]["Rate rules"])
            for family in families_sorted
        ]
    )
    ax.bar(xs, ys, bottom=bottom, label="Rate rules", edgecolor="black")
    bottom += ys

    ys = np.array(
        [
            len(rxn_by_kinetic_source_and_family[family]["QM"])
            for family in families_sorted
        ]
    )
    ax.bar(xs, ys, bottom=bottom, label="QM", edgecolor="black")
    bottom += ys

    ys = np.array(
        [
            len(rxn_by_kinetic_source_and_family[family]["Library"])
            for family in families_sorted
        ]
    )
    ax.bar(xs, ys, bottom=bottom, label="Library/Training", edgecolor="black")

    ax.set_xticks(xs)
    ax.set_xticklabels(families_sorted, rotation=45, ha="right")


(model_name,) = parse_arguments()

print("model_name: ", model_name)

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

fragments = set(["CD", "CDB", "AR", "KR", "AH", "PR", "CP", "OR", "OH"])

print("Loading mechanism files...")

if model_name in ["basecase_debutanizer_model", "QCMD_cell_model"]:
    debutanizer_liq_chemkin_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/liquid_mechanism/chem_annotated.inp"
    debutanizer_liq_spc_dict_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/liquid_mechanism/species_dictionary.txt"
    debutanizer_film_chemkin_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/film_mechanism/chem_annotated_film.inp"
    debutanizer_film_spc_dict_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/film_mechanism/species_dictionary_film.txt"

    QCMD_liq_chemkin_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/QCMD_cell_model/liquid_mechanism/chem_annotated.inp"
    QCMD_liq_spc_dict_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/QCMD_cell_model/liquid_mechanism/species_dictionary.txt"
    QCMD_film_chemkin_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/QCMD_cell_model/film_mechanism/chem_annotated_film.inp"
    QCMD_film_spc_dict_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/QCMD_cell_model/film_mechanism/species_dictionary_film.txt"

    debutanizer_liqspcs, debutanizer_liqrxns = load_chemkin_file(
        debutanizer_liq_chemkin_path, debutanizer_liq_spc_dict_path
    )
    print(
        f"Debutanizer liquid mechanism: {len(debutanizer_liqspcs)} species, {len(debutanizer_liqrxns)} reactions"
    )
    debutanizer_filmspcs, debutanizer_filmrxns = load_chemkin_file(
        debutanizer_film_chemkin_path, debutanizer_film_spc_dict_path
    )
    debutanizer_filmspcs = [
        spc
        for spc in debutanizer_filmspcs
        if spc.label in fragments or "(L)" in spc.label
    ]
    print(
        f"Debutanizer film mechanism: {len(debutanizer_filmspcs)} species, {len(debutanizer_filmrxns)} reactions"
    )

    QCMD_liqspcs, QCMD_liqrxns = load_chemkin_file(
        QCMD_liq_chemkin_path, QCMD_liq_spc_dict_path
    )
    print(
        f"QCMD liquid mechanism: {len(QCMD_liqspcs)} species, {len(QCMD_liqrxns)} reactions"
    )
    QCMD_filmspcs, QCMD_filmrxns = load_chemkin_file(
        QCMD_film_chemkin_path, QCMD_film_spc_dict_path
    )
    QCMD_filmspcs = [
        spc for spc in QCMD_filmspcs if spc.label in fragments or "(L)" in spc.label
    ]
    print(
        f"QCMD film mechanism: {len(QCMD_filmspcs)} species, {len(QCMD_filmrxns)} reactions"
    )

    all_debutanizer_spcs = debutanizer_liqspcs + debutanizer_filmspcs
    all_QCMD_spcs = QCMD_liqspcs + QCMD_filmspcs
    all_debutanizer_rxns = debutanizer_liqrxns + debutanizer_filmrxns
    all_QCMD_rxns = QCMD_liqrxns + QCMD_filmrxns
elif model_name == "trace_oxygen_perturbed_debutanizer_model":
    debutanizer_liq_chemkin_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/trace_oxygen_perturbed/liquid_mechanism/chem_annotated.inp"
    debutanizer_liq_spc_dict_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/trace_oxygen_perturbed/liquid_mechanism/species_dictionary.txt"
    debutanizer_film_chemkin_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/trace_oxygen_perturbed/film_mechanism/chem_annotated_film.inp"
    debutanizer_film_spc_dict_path = "/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/trace_oxygen_perturbed/film_mechanism/species_dictionary_film.txt"

    debutanizer_liqspcs, debutanizer_liqrxns = load_chemkin_file(
        debutanizer_liq_chemkin_path, debutanizer_liq_spc_dict_path
    )
    print(
        f"Debutanizer liquid mechanism: {len(debutanizer_liqspcs)} species, {len(debutanizer_liqrxns)} reactions"
    )
    debutanizer_filmspcs, debutanizer_filmrxns = load_chemkin_file(
        debutanizer_film_chemkin_path, debutanizer_film_spc_dict_path
    )
    debutanizer_filmspcs = [
        spc
        for spc in debutanizer_filmspcs
        if spc.label in fragments or "(L)" in spc.label
    ]
    print(
        f"Debutanizer film mechanism: {len(debutanizer_filmspcs)} species, {len(debutanizer_filmrxns)} reactions"
    )

    all_debutanizer_spcs = debutanizer_liqspcs + debutanizer_filmspcs
    all_debutanizer_rxns = debutanizer_liqrxns + debutanizer_filmrxns

print("Plotting model parameter source distributions...")
print(f"Model: {model_name}")

if model_name in ["basecase_debutanizer_model", "QCMD_cell_model"]:
    fig, axs = plt.subplots(3, 2, figsize=(8, 11))

    ax = axs[0, 0]
    spc_by_thermo_source = categorize_spc_by_thermo_source(all_debutanizer_spcs)
    plot_spc_thermo_source_by_spc_size(spc_by_thermo_source, ax=ax)
    ax.legend(bbox_to_anchor=(0.0, -0.1, 1.0, -0.1), loc=1)
    ax.set_title("(a) Debutanizer model", loc="left")
    ax.set_ylabel("Number of species")
    ax.set_xlabel("Number of heavy atoms")

    ax = axs[1, 0]
    categorize_rxn_by_kinetic_source_and_size(all_debutanizer_rxns, ax=ax)
    ax.set_title("(b) Debutanizer model", loc="left")
    ax.set_ylabel("Number of reactions")
    ax.set_xlabel("Number of heavy atoms")

    for rxn in all_debutanizer_rxns:
        rxn.family = get_reaction_family(rmg, rxn)

    ax = axs[2, 0]
    categorize_rxn_by_kinetic_source_and_family(all_debutanizer_rxns, ax=ax)
    ax.set_title("(c) Debutanizer model", loc="left")
    ax.set_ylabel("Number of reactions")
    ax.set_xlabel("Reaction family")

    ax = axs[0, 1]
    categorize_spc_by_thermo_source(all_QCMD_spcs, ax=ax)
    ax.set_title("(c) QCMD cell model", loc="left")
    ax.set_xlabel("Number of heavy atoms")

    ax = axs[1, 1]
    categorize_rxn_by_kinetic_source_and_size(all_QCMD_rxns, ax=ax)
    ax.set_title("(d) QCMD cell model", loc="left")
    ax.set_xlabel("Number of heavy atoms")

    for rxn in all_QCMD_rxns:
        rxn.family = get_reaction_family(rmg, rxn)

    ax = axs[2, 1]
    categorize_rxn_by_kinetic_source_and_family(all_QCMD_rxns, ax=ax)
    ax.set_title("(e) QCMD cell model", loc="left")
    ax.set_xlabel("Reaction family")

    fig.tight_layout()
    fig.savefig("Figures/thermo_and_kinetic_source.pdf", bbox_inches="tight")

elif model_name == "trace_oxygen_perturbed_debutanizer_model":
    fig, axs = plt.subplots(2, 3, figsize=(11, 8))
    ax = axs[0, 0]
    spc_by_thermo_source = categorize_spc_by_thermo_source(all_debutanizer_spcs)
    plot_spc_thermo_source_by_spc_size(spc_by_thermo_source, ax=ax)
    ax.legend(bbox_to_anchor=(0.0, -0.1, 1.0, -0.1), loc=1)
    ax.set_title("(a)", loc="left")
    ax.set_ylabel("Number of species")
    ax.set_xlabel("Number of heavy atoms")

    ax = axs[0, 1]
    plot_spc_thermo_source_by_spc_size(spc_by_thermo_source, ax=ax, use_carbon=True)
    ax.set_title("(b)", loc="left")
    ax.set_xlabel("Number of carbon atoms")

    ax = axs[0, 2]
    plot_spc_thermo_source_by_spc_size(spc_by_thermo_source, ax=ax, use_oxygen=True)
    ax.set_title("(c)", loc="left")
    ax.set_xlabel("Number of oxygen atoms")

    ax = axs[1, 0]
    categorize_rxn_by_kinetic_source_and_size(all_debutanizer_rxns, ax=ax)
    ax.set_title("(d)", loc="left")
    ax.set_ylabel("Number of reactions")
    ax.set_xlabel("Number of heavy atoms")

    for rxn in all_debutanizer_rxns:
        rxn.family = get_reaction_family(rmg, rxn)

    ax = axs[1, 1]
    categorize_rxn_by_kinetic_source_and_family(all_debutanizer_rxns, ax=ax)
    ax.set_title("(e)", loc="left")
    ax.set_xlabel("Reaction family")

    ax = axs[1, 2]
    ax.axis("off")

    fig.align_labels()
    fig.tight_layout()
    fig.savefig("Figures/thermo_and_kinetic_source.pdf", bbox_inches="tight")

print("Plotting model size distributions...")
if model_name in ["basecase_debutanizer_model", "QCMD_cell_model"]:
    fig, axs = plt.subplots(2, 2, figsize=(8, 8))

    ax = axs[0, 0]
    plot_num_rxns_per_family(debutanizer_liqrxns, ax=ax)
    ax.set_title("Debutanizer liquid submodel")
    ax.set_ylabel("Number of reactions")

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

elif model_name == "trace_oxygen_perturbed_debutanizer_model":
    fig, axs = plt.subplots(2, 2, figsize=(8, 7))

    ax = axs[0, 0]
    plot_num_rxns_per_family(debutanizer_liqrxns, ax=ax)
    ax.set_title("Debutanizer liquid submodel")
    ax.set_ylabel("Number of reactions")
    ax.set_xlabel("Reaction family")

    ax = axs[0, 1]
    plot_num_rxns_per_family(debutanizer_filmrxns, ax=ax)
    ax.set_title("Debutanizer film submodel")
    ax.set_xlabel("Reaction family")

    ax = axs[1, 0]
    plot_num_rxns_per_spc(debutanizer_liqspcs, debutanizer_liqrxns, ax=ax)
    ax.set_ylabel("Number of species")
    ax.set_xlabel("Number of involved reactions")

    ax = axs[1, 1]
    plot_num_rxns_per_spc(debutanizer_filmspcs, debutanizer_filmrxns, ax=ax)
    ax.set_xlabel("Number of involved reactions")

    fig.tight_layout()
    fig.savefig("Figures/num_rxns_per_family_and_per_spc.pdf", bbox_inches="tight")
