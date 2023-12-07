# %%
import os
import argparse

import yaml
import pandas as pd
from tqdm import tqdm

import numpy as np

import matplotlib.pyplot as plt

# change default font size to 12
plt.rcParams.update({"font.size": 12})

# %%
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model_name",
        type=str,
        required=True,
        help="The name of the model.",
    )
    parser.add_argument(
        "--all_simulation_directory",
        type=str,
        required=True,
        help="The path to all simulation results.",
    )

    args = parser.parse_args(
        [
            "--model_name", "trace_oxygen_perturbed", "--all_simulation_directory", "/home/gridsan/hwpang/Jobs/PolymerFoulingModeling/TraceO2Perturbed/trace_oxygen_perturbed_20230718/simulation_results"
        ]
    )
    model_name = args.model_name
    all_simulation_directory = args.all_simulation_directory

    return (
        model_name,
        all_simulation_directory,
    )


(
    model_name,
    all_simulation_directory,
) = parse_arguments()

# %%
perturbed_species_list = ["OXYGEN"]
perturbed_factor_string_list = [
    "0.0",
    "1e-6",
    "1e-5",
    "1e-4",
    "1e-3",
    "1e-2",
    "1e-1",
    "1e0",
]
perturbed_factor_num_list = [float(perturbed_factor) for perturbed_factor in perturbed_factor_string_list]
delta_t = 32.0
trays = range(1, 41)

# %%
print("Loading all alpha and rates...")
asf_results = dict()

for perturbed_species in perturbed_species_list:
    for perturbed_factor in tqdm(perturbed_factor_string_list):
        
        alpha_rates_path = os.path.join(
            all_simulation_directory,
            f"{perturbed_species}_{perturbed_factor}_3600.0_{delta_t}",
            "alpha_rates.yml"
        )
        
        
        with open(alpha_rates_path, "r") as f:
            results = yaml.load(f, Loader=yaml.FullLoader)
        all_alphas, consumption_rates, production_rates, rxn_rates, radical_labels = results
        alpha1, alpha2, alphas, alphas_DA = all_alphas
        # R_labels, ROO_labels, RO_labels = radical_labels
        alpha1 = np.array(alpha1)
        alpha2 = np.array(alpha2)
        alphas = np.array(alphas)
        # alphas_DA = np.array(alphas_DA)
        for key in consumption_rates.keys():
            consumption_rates[key] = np.array(consumption_rates[key])
        alpha1_RC_add = alpha1 * (
            consumption_rates["R._Add"]
            / (consumption_rates["R._Add"] + consumption_rates["R.+O2"])
        )
        alpha1_RC_O2 = alpha1 * (
            consumption_rates["R.+O2"]
            / (consumption_rates["R._Add"] + consumption_rates["R.+O2"])
        )
        asf_results[perturbed_species, perturbed_factor] = alpha1_RC_add, alpha1_RC_O2, alphas


# %%
factorcmap = plt.get_cmap("PuRd")

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8, 6), sharex=True, sharey=True)

for perturbed_species in perturbed_species_list:
    for perturbed_factor_ind, perturbed_factor in enumerate(perturbed_factor_string_list):
        
        last_one = perturbed_factor_ind == len(perturbed_factor_string_list)-1
        
        alpha1_RC_add, alpha1_RC_O2, alphas = asf_results[perturbed_species, perturbed_factor]
        
        ax = axs.flat[0]
        ax.scatter(
            trays,
            alpha1_RC_add,
            color=factorcmap(
                perturbed_factor_ind / len(perturbed_factor_string_list)
            ),
            zorder=2,
            marker="s",
        )
        ax.plot(trays, alpha1_RC_add, zorder=1, color="grey")
        ax.set_title("(a.1) $\\alpha_\mathrm{RC. add}$", loc="left")

        ax = axs.flat[1]
        ax.scatter(
            trays,
            alpha1_RC_O2,
            color=factorcmap(
                perturbed_factor_ind / len(perturbed_factor_string_list)
            ),
            zorder=2,
            marker="x",
        )
        ax.plot(trays, alpha1_RC_O2, zorder=1, color="grey")
        ax.set_ylim([0, 1])
        ax.set_title("(a.2) $\\alpha_\mathrm{RC.+O_2}$", loc="left")

        ax = axs.flat[2]
        ax.scatter(trays, alpha2, color=factorcmap(
                perturbed_factor_ind / len(perturbed_factor_string_list)
            ), zorder=2)
        ax.plot(trays, alpha2, zorder=1, color="grey")
        ax.set_ylim([0, 1])
        ax.set_title("(b) $\\alpha_\mathrm{ROO.}$", loc="left")

        ax = axs.flat[3]
        ax.scatter(trays, alphas, color=factorcmap(
                perturbed_factor_ind / len(perturbed_factor_string_list)
            ), zorder=2)
        ax.plot(trays, alphas, zorder=1, color="grey")
        ax.set_title("(c) $\\alpha_\mathrm{R.}$", loc="left")

fig.supylabel("(-)")
fig.supxlabel("Tray", y=0.04)

cbar_ax = fig.add_axes([1.0, 0.15, 0.02, 0.7])
sm = plt.cm.ScalarMappable(
    cmap="RdPu", norm=plt.matplotlib.colors.LogNorm(vmin=1e-7, vmax=1e0)
)
cbar = fig.colorbar(
    sm,
    ticks=perturbed_factor_num_list,
    orientation="vertical",
    label="$[\mathrm{O}_2]/[\mathrm{O}_2]_\mathrm{sat}$",
    cax=cbar_ax,
)
fig.tight_layout()
fig.savefig(f"Figures/{model_name}_ASF_distribution_sens.pdf", bbox_inches="tight")

# %%



