import os
import yaml
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils import get_rops

# change default font size to 12
plt.rcParams.update({"font.size": 12})


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

    args = parser.parse_args()
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

print("model_name: ", model_name)
print("all_simulation_directory: ", all_simulation_directory)

d = 2.5
h = 0.3
A = (d / 2) ** 2 * np.pi
Vliq = A * h

trays = range(1, 41)
if model_name == "basecase_debutanizer_model":
    perturbed_species_list = ["1,3-BUTADIENE", "CYCLOPENTADIENE"]
    perturbed_factor_list = ["0.5", "0.7", "0.9", "1.1", "1.3", "1.5", "1.7", "1.9"]
    delta_t = 64.0
    keys = ["Diels-Alder addition", "radical addition"]
else:
    perturbed_species_list = ["OXYGEN"]
    perturbed_factor_list = ["0.0", "1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1", "1e0"]
    delta_t = 32.0
    keys = ["oxygen-center radical addition", "carbon-center radical addition"]

perturbed_factor_num_list = [
    float(perturbed_factor) for perturbed_factor in perturbed_factor_list
]
tf_liq = np.arange(0, 3600.0 + delta_t, delta_t)[-1]

print("Loading all liquid simulation results...")
liquid_simulations = dict()
for perturbed_species in perturbed_species_list:
    for perturbed_factor in perturbed_factor_list:
        simulation_path = os.path.join(
            all_simulation_directory,
            f"{perturbed_species}_{perturbed_factor}_3600.0_{delta_t}",
            f"simulation_vapor_liquid_yliqn_{tf_liq}.csv",
        )
        liquid_simulations[perturbed_species, perturbed_factor] = pd.read_csv(
            simulation_path
        )

print("Loading all film simulation results...")

film_simulations = dict()
for perturbed_species in perturbed_species_list:
    for perturbed_factor in perturbed_factor_list:
        for tray in trays:
            simulation_path = os.path.join(
                all_simulation_directory,
                f"{perturbed_species}_{perturbed_factor}_3600.0_{delta_t}",
                f"simulation_film_{tray}.csv",
            )
            film_simulations[perturbed_species, perturbed_factor, tray] = pd.read_csv(
                simulation_path
            )

print("Loading all film rop results...")
film_rop_results = dict()
for perturbed_species in perturbed_species_list:
    for perturbed_factor in perturbed_factor_list:
        for tray in trays:
            film_rop_results[(perturbed_species, perturbed_factor, tray)] = pd.read_csv(
                os.path.join(
                    all_simulation_directory,
                    f"{perturbed_species}_{perturbed_factor}_3600.0_{delta_t}",
                    f"simulation_film_rop_{tray}.csv",
                )
            )

print("Loading radical labels...")

alpha_rates_path = os.path.join(
    all_simulation_directory,
    f"{perturbed_species_list[0]}_{perturbed_factor_list[0]}_3600.0_{delta_t}",
    "alpha_rates.yml",
)

with open(alpha_rates_path, "r") as f:
    results = yaml.load(f, Loader=yaml.FullLoader)
    all_alphas, consumption_rates, production_rates, rxn_rates, radical_labels = results
    alpha1, alpha2, alphas, alphas_DA = all_alphas
    R_labels, ROO_labels, RO_labels = radical_labels

# define helper functions


def calculate_film_growth_time_constant(df):
    n_rows = len(df.index)
    m = df.loc[n_rows - 2, "mass"]
    dmdt = (df.loc[n_rows - 1, "mass"] - df.loc[n_rows - 3, "mass"]) / (
        df.loc[n_rows - 1, "timestamp"] - df.loc[n_rows - 3, "timestamp"]
    )
    return m / dmdt / 3600.0 / 24.0 / 365.0  # year


def calculate_fragment_per_mass(df, label):
    n_rows = len(df.index)
    return df.loc[n_rows - 1, label] / df.loc[n_rows - 1, "mass"]


def calculate_fouling_chemistry_contribution(film_rop_results, normalize=True):

    fouling_chemistry_contribution = dict()
    for perturbed_species in perturbed_species_list:
        for perturbed_factor in perturbed_factor_list:

            fouling_chemistry_contribution[perturbed_species, perturbed_factor] = {key: np.zeros(len(trays)) for key in keys}
            
            for tray in trays:
                df = film_rop_results[perturbed_species, perturbed_factor, tray]
                name_inds = df["rop_spcname"] == "mass"
                rop_rxnstrs = df.loc[name_inds, "rop_rxnstr"]
                rop_rxncomments = df.loc[name_inds, "rop_rxncomment"]
                rops = df.loc[name_inds, "rop"]

                if model_name == "basecase_debutanizer_model":
                    for rop_rxncomment, rop in zip(rop_rxncomments, rops):
                        for key in keys:
                            if key in rop_rxncomment:
                                fouling_chemistry_contribution[perturbed_species, perturbed_factor][key][tray-1] += rop
                else:
                    for rop_rxnstr, rop_rxncomment, rop in zip(rop_rxnstrs, rop_rxncomments, rops):
                        if ("O" in rop_rxnstr and "radical addition" in rop_rxncomment) or "[O][O](L) + AR <=> PR" in rop_rxncomment or "[O][O](L) + KR <=> PR" in rop_rxncomment:
                            fouling_chemistry_contribution[perturbed_species, perturbed_factor]["oxygen-center radical addition"][tray-1] += rop
                        elif "radical addition" in rop_rxncomment:
                            fouling_chemistry_contribution[perturbed_species, perturbed_factor]["carbon-center radical addition"][tray-1] += rop
    
    if normalize:
        for perturbed_species in perturbed_species_list:
            for perturbed_factor in perturbed_factor_list:
                rop_sum = np.sum(fouling_chemistry_contribution[perturbed_species, perturbed_factor][key] for key in keys)
                for key in keys:
                    fouling_chemistry_contribution[perturbed_species, perturbed_factor][key] /= rop_sum
    
    return fouling_chemistry_contribution


print("Plotting all sensitivity simulation results...")

factorcmap = plt.get_cmap("PuRd")
nrows = 5
ncols = 2
markers = ["x", "v", "D", "P", "X", "8", "p", "*", "h", "H", "d"]
# cmaps = ["Blues", "Greens", "Purples", "Oranges"]
fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(6, 10), sharex=True)
normalized_fouling_chemistry_contribution = calculate_fouling_chemistry_contribution(film_rop_results, normalize=True)

if model_name == "basecase_debutanizer_model":
    
    for species_ind, perturbed_species in enumerate(perturbed_species_list):
        for factor_ind, perturbed_factor in enumerate(perturbed_factor_list):
            ax = axs[0, species_ind]
            liquid_carbon_center_radical_concs = (
                liquid_simulations[perturbed_species, perturbed_factor]
                .loc[:, R_labels]
                .sum(axis=1)
                / Vliq
            )
            ax.scatter(
                trays,
                liquid_carbon_center_radical_concs,
                color=factorcmap(factor_ind / len(perturbed_factor_list)),
            )
            ax.set_ylabel("R.(L) (mol/m³)")
            ax.set_yscale("log")
            ax.set_title(perturbed_species)

            ax = axs[1, species_ind]
            film_growth_time_constants = [
                calculate_film_growth_time_constant(
                    film_simulations[perturbed_species, perturbed_factor, tray]
                )
                for tray in trays
            ]
            ax.scatter(
                trays,
                film_growth_time_constants,
                color=factorcmap(factor_ind / len(perturbed_factor_list)),
            )
            ax.set_ylabel("Film growth τ (year)")
            ax.set_yscale("log")

            ax = axs[2, species_ind]
            film_fragment_per_mass = [
                calculate_fragment_per_mass(
                    film_simulations[perturbed_species, perturbed_factor, tray],
                    "AR",
                )
                for tray in trays
            ]
            ax.scatter(
                trays,
                film_fragment_per_mass,
                color=factorcmap(factor_ind / len(perturbed_factor_list)),
            )
            ax.set_ylabel("AR (mol/kg)")
            ax.set_yscale("log")

            ax = axs[3, species_ind]
            film_fragment_per_mass = [
                calculate_fragment_per_mass(
                    film_simulations[perturbed_species, perturbed_factor, tray],
                    "KR",
                )
                for tray in trays
            ]
            ax.scatter(
                trays,
                film_fragment_per_mass,
                color=factorcmap(factor_ind / len(perturbed_factor_list)),
            )
            ax.set_ylabel("KR (mol/kg)")
            ax.set_yscale("log")

            ax = axs[4, species_ind]
            contribution = normalized_fouling_chemistry_contribution[perturbed_species, perturbed_factor]
            label = None
            for ind, key in enumerate(keys):
                if perturbed_factor == perturbed_factor_list[-1] and perturbed_species == perturbed_species_list[-1]:
                    label = key
                ax.scatter(
                    trays,
                    contribution[key],
                    # color=plt.get_cmap(cmaps[ind])(factor_ind / len(perturbed_factor_list)),
                    color=factorcmap(factor_ind / len(perturbed_factor_list)),
                    marker=markers[ind],
                    edgecolors="black",
                    label=label,
                    zorder=ind,
                )
            ax.set_ylabel("Film growth\nchemistry (%)")

            if perturbed_factor == perturbed_factor_list[-1] and perturbed_species == perturbed_species_list[-1]:
                ax.legend(bbox_to_anchor=(1.05, -0.3), loc="upper right")

else:
    for species_ind, perturbed_species in enumerate(perturbed_species_list):
        for factor_ind, perturbed_factor in enumerate(perturbed_factor_list):
            ax = axs[0, 0]
            liquid_carbon_center_radical_concs = (
                liquid_simulations[perturbed_species, perturbed_factor]
                .loc[:, R_labels]
                .sum(axis=1)
                / Vliq
            )
            ax.scatter(
                trays,
                liquid_carbon_center_radical_concs,
                color=factorcmap(factor_ind / len(perturbed_factor_list)),
            )
            ax.set_ylabel("RC.(L) (mol/m³)")
            ax.set_yscale("log")
            ax.set_title("(a)", loc="left")

            if perturbed_factor != "0.0":
                ax = axs[0, 1]
                liquid_peroxyl_radical_concs = (
                    liquid_simulations[perturbed_species, perturbed_factor]
                    .loc[:, ROO_labels]
                    .sum(axis=1)
                    / Vliq
                )
                ax.scatter(
                    trays,
                    liquid_peroxyl_radical_concs,
                    color=factorcmap(factor_ind / len(perturbed_factor_list)),
                )
                ax.set_ylabel("ROO.(L) (mol/m³)")
                ax.set_yscale("log")
                ax.set_title("(b)", loc="left")

                ax = axs[1, 0]
                liquid_alkoxyl_radical_concs = (
                    liquid_simulations[perturbed_species, perturbed_factor]
                    .loc[:, RO_labels]
                    .sum(axis=1)
                    / Vliq
                )
                ax.scatter(
                    trays,
                    liquid_alkoxyl_radical_concs,
                    color=factorcmap(factor_ind / len(perturbed_factor_list)),
                )
                ax.set_ylabel("RO.(L) (mol/m³)")
                ax.set_yscale("log")
                ax.set_title("(c)", loc="left")

            ax = axs[1, 1]
            film_growth_time_constants = [
                calculate_film_growth_time_constant(
                    film_simulations[perturbed_species, perturbed_factor, tray]
                )
                for tray in trays
            ]
            ax.scatter(
                trays,
                film_growth_time_constants,
                color=factorcmap(factor_ind / len(perturbed_factor_list)),
            )
            ax.set_ylabel("Film growth τ (yr)")
            ax.set_yscale("log")
            ax.set_title("(d)", loc="left")

            ax = axs[2, 0]
            AR_concs = [
                calculate_fragment_per_mass(
                    film_simulations[perturbed_species, perturbed_factor, tray], "AR"
                )
                for tray in trays
            ]
            ax.scatter(
                trays, AR_concs, color=factorcmap(factor_ind / len(perturbed_factor_list))
            )
            ax.set_ylabel("AR/mass (mol/kg)")
            ax.set_yscale("log")
            ax.set_title("(e)", loc="left")

            ax = axs[2, 1]
            KR_concs = [
                calculate_fragment_per_mass(
                    film_simulations[perturbed_species, perturbed_factor, tray], "KR"
                )
                for tray in trays
            ]
            ax.scatter(
                trays, KR_concs, color=factorcmap(factor_ind / len(perturbed_factor_list))
            )
            ax.set_ylabel("KR/mass (mol/kg)")
            ax.set_yscale("log")
            ax.set_title("(f)", loc="left")

            if perturbed_factor != "0.0":
                ax = axs[3, 0]
                PR_concs = [
                    calculate_fragment_per_mass(
                        film_simulations[perturbed_species, perturbed_factor, tray], "PR"
                    )
                    for tray in trays
                ]
                ax.scatter(
                    trays,
                    PR_concs,
                    color=factorcmap(factor_ind / len(perturbed_factor_list)),
                )
                ax.set_ylabel("PR/mass (mol/kg)")
                ax.set_yscale("log")
                ax.set_title("(g)", loc="left")

                ax = axs[3, 1]
                OR_concs = [
                    calculate_fragment_per_mass(
                        film_simulations[perturbed_species, perturbed_factor, tray], "OR"
                    )
                    for tray in trays
                ]
                ax.scatter(
                    trays,
                    OR_concs,
                    color=factorcmap(factor_ind / len(perturbed_factor_list)),
                )
                ax.set_ylabel("OR/mass (mol/kg)")
                ax.set_yscale("log")
                ax.set_title("(h)", loc="left")

                ax = axs[4, 0]
                contribution = normalized_fouling_chemistry_contribution[perturbed_species, perturbed_factor]
                label = None
                for ind, key in enumerate(keys):
                    if perturbed_factor == perturbed_factor_list[-1] and perturbed_species == perturbed_species_list[-1]:
                        label = key
                    ax.scatter(
                        trays,
                        contribution[key],
                        color=factorcmap(factor_ind / len(perturbed_factor_list)),
                        marker=markers[ind],
                        edgecolors="black",
                        label=label,
                        zorder=ind,
                    )
                ax.set_ylabel("Film growth\nchemistry (%)")

                if perturbed_factor == perturbed_factor_list[-1] and perturbed_species == perturbed_species_list[-1]:
                    ax.legend(bbox_to_anchor=(1.05, -0.3), loc="upper right")

axs[-1, 0].set_xlabel("Trays")
axs[-1, 1].set_xlabel("Trays")

cbar_ax = fig.add_axes([1.0, 0.15, 0.02, 0.7])
if model_name == "basecase_debutanizer_model":
    sm = plt.cm.ScalarMappable(
        cmap="RdPu", norm=plt.Normalize(vmin=0.5, vmax=1.9)
    )
else:
    sm = plt.cm.ScalarMappable(
        cmap="RdPu", norm=plt.matplotlib.colors.LogNorm(vmin=1e-7, vmax=1e0)
    )
cbar = fig.colorbar(
    sm,
    ticks=perturbed_factor_num_list,
    orientation="vertical",
    label="Perturbation",
    cax=cbar_ax,
)

fig.tight_layout()
fig.savefig(f"Figures/{model_name}_sens_all.pdf", bbox_inches="tight")
plt.close()