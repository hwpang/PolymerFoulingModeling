import os
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
import argparse

from rmgpy import settings
from rmgpy.species import Species
from rmgpy.data.thermo import ThermoDatabase


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

    return model_name


model_name = parse_arguments()

# load thermo libraries containing QM calculated data
thermo_db = ThermoDatabase()
thermo_db.load(
    path=os.path.join(settings["database.directory"], "thermo"),
    libraries=[
        "hwpang_fouling",
        # 'thermo_combined_new',
        "Conjugated_diene",
        # 'ExptModelTraceO2',
        # 'multi_trays_combined',
        # 'multi_trays_v3',
        # 'multi_trays',
        # 'thermo_combined',
        # 's3_5_7_ane',
        # 'Xiaorui_thermo',
    ],
    depository=False,
)

QM_spc_smiles_set = set()
QM_spcs = []
smis = []

for library_label, library in thermo_db.libraries.items():
    print(library_label, len(library.entries))
    count = 0
    for entry in library.entries.values():
        mol = entry.item
        spc = Species(molecule=[mol])
        spc.generate_resonance_structures()
        spc.thermo = entry.data
        smi = spc.molecule[0].to_smiles()
        if model_name == "basecase_debutanizer_model":
            if "N" in smi or "n" in smi:
                continue
        elif model_name == "trace_oxygen_perturbed_debutanizer_model":
            if "N" in smi or "n" in smi:
                continue
            if "O" not in smi and "o" not in smi:
                continue
        else:
            raise ValueError("Unknown model name.")

        if smi not in QM_spc_smiles_set:
            QM_spc_smiles_set.add(smi)
            QM_spcs.append(spc)
            if library_label != "thermo_combined_new":
                smis.append(smi)
            count += 1
    print(count)

print(smis)
print(len(smis))

# esitmate thermo with GAV
GAV_spcs = []
for spc in QM_spcs:
    new_spc = deepcopy(spc)
    new_spc.thermo = thermo_db.estimate_thermo_via_group_additivity(new_spc.molecule[0])
    GAV_spcs.append(new_spc)

os.makedirs("Figures", exist_ok=True)

# plot comparison
QM_H298 = np.array([spc.thermo.get_enthalpy(298) for spc in QM_spcs]) / 1000 / 4.184
GAV_H298 = np.array([spc.thermo.get_enthalpy(298) for spc in GAV_spcs]) / 1000 / 4.184
QM_S298 = np.array([spc.thermo.get_entropy(298) for spc in QM_spcs]) / 4.184
GAV_S298 = np.array([spc.thermo.get_entropy(298) for spc in GAV_spcs]) / 4.184

# categorize the species by number of carbon atoms and oxygen atoms
spc_ind_by_carbon_num = {}
spc_ind_by_oxygen_num = {}

for ind, spc in enumerate(QM_spcs):
    num = spc.molecule[0].get_num_atoms("C")
    if num == 0:
        continue
    if num not in spc_ind_by_carbon_num:
        spc_ind_by_carbon_num[num] = []
    spc_ind_by_carbon_num[num].append(ind)

    num = spc.molecule[0].get_num_atoms("O")
    if num not in spc_ind_by_oxygen_num:
        spc_ind_by_oxygen_num[num] = []
    spc_ind_by_oxygen_num[num].append(ind)

H298_diff_by_carbon_num = {}
S298_diff_by_carbon_num = {}
for num, ind_list in spc_ind_by_carbon_num.items():
    H298_diff_by_carbon_num[num] = QM_H298[ind_list] - GAV_H298[ind_list]
    S298_diff_by_carbon_num[num] = QM_S298[ind_list] - GAV_S298[ind_list]

H298_diff_by_oxygen_num = {}
S298_diff_by_oxygen_num = {}
for num, ind_list in spc_ind_by_oxygen_num.items():
    H298_diff_by_oxygen_num[num] = QM_H298[ind_list] - GAV_H298[ind_list]
    S298_diff_by_oxygen_num[num] = QM_S298[ind_list] - GAV_S298[ind_list]

# plot violin plot
min_carbon_num = min(H298_diff_by_carbon_num.keys())
max_carbon_num = max(H298_diff_by_carbon_num.keys())
min_oxygen_num = min(H298_diff_by_oxygen_num.keys())
max_oxygen_num = max(H298_diff_by_oxygen_num.keys())

if model_name == "basecase_debutanizer_model":
    fig, axes = plt.subplots(1, 2, figsize=(8, 3))
    ax = axes[0]
    ax.hist(QM_H298 - GAV_H298, edgecolor="black", bins=50)
    ax.set_xlabel("$\Delta_\mathrm{f}H^o$(298 K) QM - GAV (kcal/mol)")
    ax.set_ylabel("Count")

    ax = axes[1]
    ax.hist(QM_S298 - GAV_S298, edgecolor="black", bins=50)
    ax.set_xlabel("$\Delta_\mathrm{f}S^o$(298 K) QM - GAV (cal/mol/K)")
    ax.set_ylabel("Count")

    fig.tight_layout()
    fig.savefig("Figures/QM_GAV_comparison.pdf")

    fig, axes = plt.subplots(2, 2, figsize=(8, 6))

    ax = axes[0, 0]
    ax.violinplot(
        list(H298_diff_by_carbon_num.values()),
        positions=list(H298_diff_by_carbon_num.keys()),
        showmeans=True,
    )
    ax.plot([min_carbon_num - 0.5, max_carbon_num + 0.5], [0, 0], "--", color="gray")
    ax.set_xticks(range(min_carbon_num, max_carbon_num + 1))
    ax.set_xticklabels(range(min_carbon_num, max_carbon_num + 1))
    ax.set_ylabel("$\Delta_\mathrm{f}H^o$(298 K) QM - GAV (kcal/mol)")

    ax = axes[1, 0]
    ax.violinplot(
        list(S298_diff_by_carbon_num.values()),
        positions=list(S298_diff_by_carbon_num.keys()),
        showmeans=True,
    )
    ax.plot([min_carbon_num - 0.5, max_carbon_num + 0.5], [0, 0], "--", color="gray")
    ax.set_xticks(range(min_carbon_num, max_carbon_num + 1))
    ax.set_xticklabels(range(min_carbon_num, max_carbon_num + 1))
    ax.set_ylabel("$\Delta_\mathrm{f}S^o$(298 K) QM - GAV (cal/mol/K)")
    ax.set_xlabel("Number of carbon atoms")

    ax = axes[0, 1]
    ax.violinplot(
        list(H298_diff_by_oxygen_num.values()),
        positions=list(H298_diff_by_oxygen_num.keys()),
        showmeans=True,
    )
    ax.plot([min_oxygen_num - 0.5, max_oxygen_num + 0.5], [0, 0], "--", color="gray")
    ax.set_xticks(range(min_oxygen_num, max_oxygen_num + 1))
    ax.set_xticklabels(range(min_oxygen_num, max_oxygen_num + 1))

    ax = axes[1, 1]
    ax.violinplot(
        list(S298_diff_by_oxygen_num.values()),
        positions=list(S298_diff_by_oxygen_num.keys()),
        showmeans=True,
    )
    ax.plot([min_oxygen_num - 0.5, max_oxygen_num + 0.5], [0, 0], "--", color="gray")
    ax.set_xticks(range(min_oxygen_num, max_oxygen_num + 1))
    ax.set_xticklabels(range(min_oxygen_num, max_oxygen_num + 1))
    ax.set_xlabel("Number of oxygen atoms")

    fig.tight_layout()
    fig.savefig("Figures/QM_GAV_comparison_by_carbon_and_oxygen_num.pdf")

elif model_name == "trace_oxygen_perturbed_debutanizer_model":
    fig, axes = plt.subplots(2, 2, figsize=(8, 6))

    ax = axes[0, 0]
    ax.hist(QM_H298 - GAV_H298, edgecolor="black", bins=50)
    ax.set_title("(a) $\Delta_\mathrm{f}H^o$(298 K)", loc="left")
    ax.set_xlabel("QM - GAV (kcal/mol)")
    ax.set_ylabel("Count")

    ax = axes[0, 1]
    ax.hist(QM_S298 - GAV_S298, edgecolor="black", bins=50)
    ax.set_title("(b) $\Delta_\mathrm{f}S^o$(298 K)", loc="left")
    ax.set_xlabel("QM - GAV (cal/mol/K)")
    ax.set_ylabel("Count")

    ax = axes[1, 0]
    ax.violinplot(
        list(H298_diff_by_oxygen_num.values()),
        positions=list(H298_diff_by_oxygen_num.keys()),
        showmeans=True,
    )
    ax.plot([min_oxygen_num - 0.5, max_oxygen_num + 0.5], [0, 0], "--", color="gray")
    ax.set_xticks(range(min_oxygen_num, max_oxygen_num + 1))
    ax.set_xticklabels(range(min_oxygen_num, max_oxygen_num + 1))
    ax.set_title("(a) grouped by number of oxygen atoms", loc="left")
    ax.set_ylabel("QM - GAV (kcal/mol)")
    ax.set_xlabel("Number of oxygen atoms")

    ax = axes[1, 1]
    ax.violinplot(
        list(S298_diff_by_oxygen_num.values()),
        positions=list(S298_diff_by_oxygen_num.keys()),
        showmeans=True,
    )
    ax.plot([min_oxygen_num - 0.5, max_oxygen_num + 0.5], [0, 0], "--", color="gray")
    ax.set_xticks(range(min_oxygen_num, max_oxygen_num + 1))
    ax.set_xticklabels(range(min_oxygen_num, max_oxygen_num + 1))
    ax.set_title("(b) grouped by number of oxygen atoms", loc="left")
    ax.set_ylabel("QM - GAV (cal/mol/K)")
    ax.set_xlabel("Number of oxygen atoms")

    fig.align_labels()
    fig.tight_layout()
    fig.savefig("Figures/QM_GAV_comparison.pdf", bbox_inches="tight")

else:
    raise ValueError("Unknown model name.")
